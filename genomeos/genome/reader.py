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

Openness is not enough where the chromatin says otherwise. A bivalent or
Polycomb promoter is open in a DNase assay and silent: H1 opens every HOXA
promoter and carries H3K27me3 on ten of eleven. Where the epigenome layer
(`genome/epigenome.py`) has H3K27me3 and H3K27ac peaks for the cell and the
chromosome, an open promoter with an H3K27me3 peak and no H3K27ac peak is
called poised and is not counted as read. A promoter the DNase file calls closed
but that carries H3K4me3 and H3K27ac is called read by its marks
(`read_by_marks`): a shallow DNase experiment misses active promoters.
Measured RNA backs both rules (scripts/epigenome.py reader-check,
`reader_poised_check`). Over the genome in K562, HepG2, GM12878 and IMR-90,
poised genes are expressed as rarely as closed ones. Genes read by their marks
are expressed about as often as genes read by openness, and GM12878's 1,095 of
them more often.
"""

from __future__ import annotations

import bisect
import gzip
import io
import json
import math
import re
import urllib.parse
import urllib.request
from collections.abc import Iterable
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


POISED_MARKS = ("H3K27me3", "H3K27ac", "H3K4me3")


def _mark_index(cell_type: str, mark: str, chrom: str) -> PeakIndex | None:
    """A histone mark's peaks from the epigenome layer's cache, or None when never read."""
    from genomeos.genome import epigenome

    pk = epigenome.load_mark_peaks(cell_type, mark, chrom)
    return PeakIndex(pk) if pk is not None else None


# The open-node call, defined once. genomeos/attribution/candidates.py re-derived this expression
# line for line until 2026-09-22; it now calls these, so the two cannot drift apart.
NODE_OPEN_MIN_DENSITY = 1.0  # peaks per 100 kb: the floor below which a median split is meaningless

NODE_OPEN_BASIS = (
    "a node is open where its DNase peak density is at or above the median over that biosample's "
    f"own nodes on this chromosome, floored at {NODE_OPEN_MIN_DENSITY} peaks per 100 kb. That is a "
    "rank within one cell, not an absolute call: half of every cell's nodes are open by "
    "construction, so a difference between two cells says one ranks the node above its own median "
    "and the other does not, never that the node is shut in either. Below the floor the call stops "
    "being a rank and becomes an absolute threshold on a density, which is assay depth; over the "
    "biosamples on disk that happens only on chrY. See docs/NODES-READER-WRITER.md, 'The second "
    "copy of that threshold, in the attribution candidates'."
)


def node_open_threshold(densities: Iterable[float]) -> float:
    """The peaks-per-100 kb at or above which a node counts as open for one biosample on one
    chromosome: the median over that biosample's own nodes, floored. See NODE_OPEN_BASIS for what
    the call does and does not support."""
    d = sorted(densities)
    return max(NODE_OPEN_MIN_DENSITY, d[len(d) // 2]) if d else NODE_OPEN_MIN_DENSITY


# The between-biosample open-node reading, registered 2026-09-27 before it was computed
# (docs/NODES-READER-WRITER.md, "The reader family normalised, to close milestone 1.1").
# One absolute density for every biosample: the pooled median of peaks per 100 kb over all
# 260,026 (node, biosample) calls on disk that day, frozen here so adding a biosample does not
# redefine it. The count at this density is depth; its log-depth residual is the reading.
NODE_OPEN_REFERENCE_DENSITY = 4.79

# normalised reading -> (raw reading it normalises, covariate it is normalised against). Each is
# the standardised residual of raw ~ a + b * ln(covariate) over a panel of biosamples, so it
# exists only for a panel, never for one row: see normalise_family.
NORMALISED_READINGS = {
    "nodes_open_depth_residual": ("nodes_open_at_reference", "peaks"),
    "genes_poised_mark_residual": ("genes_poised", "h3k27me3_peaks"),
    "genes_read_depth_residual": ("genes_read", "peaks"),
    "genes_read_open_depth_residual": ("genes_read_open", "peaks"),
}


def log_fit(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    """Least squares y = a + b * ln(x); returns (a, b, s) with s = sqrt(SSR / (n - 2))."""
    lx = [math.log(x) for x in xs]
    n = len(lx)
    mx, my = sum(lx) / n, sum(ys) / n
    sxx = sum((v - mx) ** 2 for v in lx)
    b = sum((u - mx) * (v - my) for u, v in zip(lx, ys, strict=True)) / sxx if sxx else 0.0
    a = my - b * mx
    ssr = sum((v - (a + b * u)) ** 2 for u, v in zip(lx, ys, strict=True))
    s = math.sqrt(ssr / (n - 2)) if n > 2 else 0.0
    # an exact fit leaves rounding error, not a residual: call it zero rather than divide by it
    return a, b, 0.0 if s <= 1e-9 * (1.0 + max(abs(v) for v in ys)) else s


def depth_residuals(xs: list[float], ys: list[float]) -> list[float]:
    """Standardised residuals of y ~ a + b * ln(x) over the panel: how many standard errors each
    biosample sits above (+) or below (-) what its covariate predicts."""
    a, b, s = log_fit(xs, ys)
    return [(y - (a + b * math.log(x))) / s if s else 0.0 for x, y in zip(xs, ys, strict=True)]


def loo_residuals(xs: list[float], ys: list[float]) -> list[float]:
    """Each biosample's standardised residual from a fit on the others. Unlike the in-sample
    residual it is not uncorrelated with ln(x) by construction, so it is the honest test of
    whether depth was removed."""
    out = []
    for i, (x, y) in enumerate(zip(xs, ys, strict=True)):
        a, b, s = log_fit(xs[:i] + xs[i + 1 :], ys[:i] + ys[i + 1 :])
        out.append((y - (a + b * math.log(x))) / s if s else 0.0)
    return out


def normalise_family(totals: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """The normalised readings for a panel of biosamples, from their genome-wide totals (each must
    carry the raw reading and its covariate, as NORMALISED_READINGS names them). A reading whose
    raw value or covariate is missing for any biosample is None for all of them: a residual is a
    property of the panel, and a panel with a hole is a different panel.

    Beside each reading, `<name>_at_edge` is True for the biosamples at the panel's highest and
    lowest covariate. Their residual is an extrapolation that rests on the fit's form: on
    2026-09-27 the log fit put K562, alone at the top of DNase depth, below HepG2 on read share,
    and cutting every biosample to the same number of peaks put it above. Compare an edge
    biosample on a depth-matched count (scripts/reader_depth_family.py --rarefied), not on this."""
    cells = list(totals)
    out: dict[str, dict[str, Any]] = {c: {} for c in cells}
    for name, (raw, cov) in NORMALISED_READINGS.items():
        ys = [totals[c].get(raw) for c in cells]
        xs = [totals[c].get(cov) for c in cells]
        ok = len(cells) > 2 and all(v is not None for v in ys) and all(v for v in xs)
        z = depth_residuals([float(v) for v in xs], [float(v) for v in ys]) if ok else [None] * len(cells)
        edge = (
            {max(cells, key=lambda c: totals[c][cov]), min(cells, key=lambda c: totals[c][cov])}
            if ok
            else set()
        )
        for c, v in zip(cells, z, strict=True):
            out[c][name] = round(v, 4) if v is not None else None
            out[c][name + "_at_edge"] = c in edge
    return out


def read_chromosome(
    cell_type: str, chrom: str, annotation, domains: list, ccres: list, marks: bool = True
) -> dict[str, Any]:
    """What this cell type reads on one chromosome: open nodes, read genes, active enhancers."""
    idx = PeakIndex(load_peaks(cell_type, chrom))
    if not idx.peaks:
        raise FileNotFoundError(f"no peaks for {cell_type} on {chrom}; run fetch_peaks first")
    genes = [g for g in annotation.genes.values() if g.locus.chrom == chrom]
    k27me3, k27ac, k4me3 = (_mark_index(cell_type, m, chrom) if marks else None for m in POISED_MARKS)
    use_marks = k27me3 is not None and k27ac is not None
    read, silent, poised, by_marks = [], [], [], []
    for g in genes:
        if g.type != "protein_coding":
            continue
        tss = g.locus.end - 1 if g.locus.strand is Strand.MINUS else g.locus.start
        lo, hi = tss - PROMOTER_WINDOW, tss + PROMOTER_WINDOW
        hits = idx.overlapping(lo, hi)
        row = {"gene": g.symbol, "signal": round(max((v for _, _, v in hits), default=0.0), 2)}
        if not hits:
            if use_marks and k4me3 is not None and k4me3.overlapping(lo, hi) and k27ac.overlapping(lo, hi):
                by_marks.append(row)
            else:
                silent.append(row)
        elif use_marks and k27me3.overlapping(lo, hi) and not k27ac.overlapping(lo, hi):
            poised.append(row)
        else:
            read.append(row)
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
    open_density = node_open_threshold(n["peaks_per_100kb"] for n in nodes)
    open_nodes = [n for n in nodes if n["peaks_per_100kb"] >= open_density]
    open_at_reference = sum(1 for n in nodes if n["peaks_per_100kb"] >= NODE_OPEN_REFERENCE_DENSITY)
    silent_nodes = [n for n in nodes if n["peaks"] == 0]
    read_open = len(read)
    read.extend(by_marks)
    coding = len(read) + len(silent) + len(poised)
    return {
        "cell_type": cell_type,
        "chrom": chrom,
        "peaks": len(idx.peaks),
        "coding_genes": coding,
        "genes_read": len(read),
        "genes_read_open": read_open + len(poised),
        "genes_read_by_marks": len(by_marks) if use_marks else None,
        "genes_poised": len(poised) if use_marks else None,
        "genes_silent": len(silent) + len(poised),
        "genes_closed": len(silent),
        "read_fraction": round(len(read) / max(1, coding), 4),
        "top_read": read[:25],
        "_read_all": [r["gene"] for r in read],
        # every gene not read, closed or poised: consumers (blocks, report, decompile) take "not in this
        # list" as read, so a truncated list or a poised gene left out of it is a wrong answer
        "silent_genes": [s["gene"] for s in silent] + sorted(p["gene"] for p in poised),
        "poised_genes": sorted(p["gene"] for p in poised),
        "read_by_marks": sorted(r["gene"] for r in by_marks),
        "marks_used": use_marks,
        # the covariate genes_poised is normalised against (the mark that calls it), carried per row
        "h3k27me3_peaks": len(k27me3.peaks) if use_marks else None,
        "enhancers": len(enh),
        "enhancers_active": active,
        "enhancers_active_fraction": round(active / max(1, len(enh)), 4),
        # enhancers_active is dominated by how deeply this biosample's DNase experiment was
        # sequenced (Spearman 0.83 against peak count over thirteen biosamples, and an
        # out-of-sample fit put testis at z = -0.11), so the rate is reported beside the count.
        # See docs/NODES-READER-WRITER.md, "Is the rest of the reader's family assay depth too?"
        "enhancers_active_per_100k_peaks": round(active / max(1, len(idx.peaks)) * 100_000, 1),
        "nodes": len(nodes),
        "nodes_open": len(open_nodes),
        # the between-biosample counterpart: one absolute density for every biosample. The count is
        # depth; normalise_family turns it into nodes_open_depth_residual over a panel
        "nodes_open_at_reference": open_at_reference,
        "nodes_silent": len(silent_nodes),
        "silent_node_ids": [n["id"] for n in silent_nodes][:100],
        "node_table": nodes,
        "evidence": {
            "peaks": EVIDENCE,
            "read": (
                "inferred: promoter (TSS ± 1 kb) overlaps a DNase peak and is not poised (an H3K27me3 peak "
                "without an H3K27ac peak), or carries H3K4me3 and H3K27ac peaks where the DNase file has "
                "none (ENCODE Histone ChIP-seq); open or marked is necessary for transcription, not proof "
                "of it"
                if use_marks
                else "inferred: promoter (TSS ± 1 kb) overlaps a DNase peak; open is necessary for "
                "transcription, not proof of it (no H3K27me3/H3K27ac peaks read for this cell and "
                "chromosome, so a poised promoter counts as read)"
            ),
            "nodes": "inferred: CTCF domains",
            # measured over thirteen biosamples, 2026-09-22 (data/results/reader_depth_family.json):
            # how much of each reading is the assay rather than the cell type
            "depth": (
                "enhancers_active and enhancers_active_fraction are depth-dominated (Spearman 0.83 "
                "against DNase peak count; use enhancers_active_per_100k_peaks for a comparison "
                "between biosamples). genes_read_by_marks is inverse DNase depth by construction: "
                "it counts promoters the marks rescue where the DNase file called nothing. "
                "nodes_open is a median split of this biosample against itself and therefore "
                "returns half the nodes for every biosample (9,988 to 10,016 of 20,002 across "
                "thirteen); it distinguishes nothing and is withdrawn as a between-biosample "
                "reading. genes_poised tracks the H3K27me3 peak call at 0.49, so a biosample with "
                "an under-called broad mark reads as un-poised. genes_read, genes_read_open and "
                "nodes_silent are banded as partly independent of depth (|rho| 0.53 to 0.59). "
                # normalised 2026-09-27 against a registration committed before the run
                # (data/results/reader_normalised.json, reader_rarefied.json)
                "Between biosamples, read the normalised counterparts, which exist for a panel "
                "and not for one row (genomeos.genome.reader.normalise_family): "
                "nodes_open_depth_residual is the log-DNase-depth residual of "
                "nodes_open_at_reference, the nodes at or above one fixed 4.79 peaks per 100 kb; "
                "genes_poised_mark_residual is genes_poised against the H3K27me3 peak count "
                "(h3k27me3_peaks) that calls it; genes_read_depth_residual and "
                "genes_read_open_depth_residual are the read counts against DNase depth. All four "
                "replicate across odd and even autosomes (Spearman 0.92 to 1.00) and sit inside "
                "|rho| < 0.4 of their covariate when each biosample is left out of its own fit. "
                "The biosamples at the top and bottom of the covariate are extrapolations "
                "(<name>_at_edge): compare those on a depth-matched count instead. nodes_open "
                "itself stays the within-biosample rank node_open_threshold defines."
            ),
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
