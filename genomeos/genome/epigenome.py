# SPDX-License-Identifier: AGPL-3.0-or-later
"""The epigenome layer: what is measured on the chromatin, per locus and per cell type.

Reader v1 (`genome/reader.py`) says whether a locus is open in a cell (ENCODE DNase-seq
peaks). Openness is necessary for function and says nothing about which function: an
open element can carry an active-enhancer mark in one cell and a Polycomb mark in
another. This module gives the chromatin its own measured layer, read from ENCODE for
the reader's eleven biosamples:

    H3K4me3   promoter
    H3K27ac   active enhancer or active promoter
    H3K4me1   enhancer, primed or active
    H3K27me3  Polycomb repression
    H3K9me3   constitutive heterochromatin
    WGBS      DNA methylation per CpG (fraction and read coverage)

Sources, one released GRCh38 experiment per (biosample, mark), chosen by a fixed rule
and pinned in `data/results/epigenome_manifest.json`:

- Histone ChIP-seq, the experiment's default (ENCODE4 pipeline) analysis: the replicated
  narrowPeak file (pseudoreplicated where no replicated one exists), streamed once and
  kept per chromosome under `data/knowledge/epigenome/peaks`, and the pooled "fold change
  over control" bigWig, read by range per chromosome into 200 bp bin means under
  `data/knowledge/epigenome/signal`. Nothing raw is kept.
- WGBS: the "methylation state at CpG" bigBed (bedMethyl: one call per strand with its
  read coverage and percent). The per-strand methylation bigWigs ENCODE also releases
  carry every cytosine, not only CpGs (positions in CCCTT and AACTT contexts on chr21),
  so a mean over them dilutes the fraction; the bigBed is the clean source. It is read
  by range with the bigWig range machinery (`attribution/bigwig.py`, same R-tree and
  chromosome tree) and distilled into 200 bp bins: CpG calls, calls with at least
  MIN_COVERAGE reads, the sum of their per-call fractions and their reads.

The layer is a record per locus per cell type (`epigenome_at`): openness, each mark
(peak call, peak signal, fold change), methylation with coverage, and a chromatin state
inferred from the marks by fixed rules, every field with its evidence and confidence and
UNKNOWN (with the reason) where the biosample has no experiment or the locus was not read.
It is a data layer; no BioIR construct reads it yet.
"""

from __future__ import annotations

import bisect
import gzip
import hashlib
import io
import json
import math
import struct
import time
import urllib.error
import urllib.parse
import urllib.request
import zlib
from array import array
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from genomeos.attribution import bigwig as _bw
from genomeos.genome import reader

ENCODE = "https://www.encodeproject.org"
CACHE = Path("data/knowledge/epigenome")
RESULTS = Path("data/results")
MANIFEST = "epigenome_manifest"
USER_AGENT = "GenomeOS/0.9 (epigenome layer; released GRCh38 files read by range)"

#: the reader's biosamples (scripts/reader_genome_wide.py DEFAULT_CELLS). The first eleven are
#: cultured -- six lines, two differentiated from a stem line, three primary cells; testis and
#: ovary (2026-09-22) are the first bulk tissue and the first gonadal lineage in the layer, and
#: carry peaks and methylation but no fold-change signal profiles yet.
CELL_TYPES = (
    "K562",
    "HepG2",
    "GM12878",
    "H1",
    "IMR-90",
    "SK-N-SH",
    "cardiac muscle cell",
    "keratinocyte",
    "hepatocyte",
    "astrocyte",
    "CD14-positive monocyte",
    "testis",
    "ovary",
)
#: biosamples that are bulk post-mortem tissue rather than a cultured population: a single donor,
#: a mixture of cell types, and only pseudoreplicated peak calls on this portal.
TISSUE_CELL_TYPES = ("testis", "ovary")
MARKS = {
    "H3K4me3": "promoter",
    "H3K27ac": "active enhancer or active promoter",
    "H3K4me1": "enhancer, primed or active",
    "H3K27me3": "Polycomb repression",
    "H3K9me3": "constitutive heterochromatin",
}
PEAK_OUTPUTS = ("replicated peaks", "pseudoreplicated peaks")
SIGNAL_OUTPUT = "fold change over control"
WGBS_OUTPUT = "methylation state at CpG"
BIN = 200
MIN_COVERAGE = 5
PROMOTER_WINDOW = reader.PROMOTER_WINDOW
BIGBED_MAGIC = 0x8789F2EB
CHROMS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
UNKNOWN = "UNKNOWN"

EVIDENCE_MARK = "experimental: ENCODE Histone ChIP-seq (GRCh38, default analysis)"
EVIDENCE_METHYLATION = "experimental: ENCODE WGBS, methylation state at CpG (GRCh38 bedMethyl)"
EVIDENCE_STATE = (
    "inferred: fixed rules over the five measured marks' peak calls (bivalent, active promoter, "
    "active enhancer, poised or primed enhancer, Polycomb, heterochromatin); not a ChromHMM model"
)


def slug(cell_type: str) -> str:
    return reader.slug(cell_type)


# ------------------------------------------------------------------------------------------
# ENCODE portal: which experiment and which files per (biosample, mark)
# ------------------------------------------------------------------------------------------


def _portal(path: str, timeout: int = 60, retries: int = 4) -> dict:
    """GET a portal JSON; an empty search is a 404 there and comes back as an empty graph.

    Every answer is cached under data/knowledge/epigenome/portal, so a manifest build that
    is interrupted resumes where it stopped instead of asking the portal again."""
    key = CACHE / "portal" / (hashlib.sha1(path.encode()).hexdigest() + ".json.gz")
    if key.exists():
        with gzip.open(key, "rt") as fh:
            return json.load(fh)
    d = _portal_get(path, timeout, retries)
    key.parent.mkdir(parents=True, exist_ok=True)
    tmp = key.with_suffix(".tmp")
    with gzip.open(tmp, "wt") as fh:
        json.dump(d, fh)
    tmp.replace(key)
    return d


def _portal_get(path: str, timeout: int, retries: int) -> dict:
    req = urllib.request.Request(
        ENCODE + path, headers={"Accept": "application/json", "User-Agent": USER_AGENT}
    )
    delay = 2.0
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {"@graph": []}
            if attempt == retries - 1:
                raise
        except (OSError, urllib.error.URLError):  # pragma: no cover - network
            if attempt == retries - 1:
                raise
        time.sleep(delay)
        delay *= 2
    return {"@graph": []}  # pragma: no cover


def _treated(summary: str) -> bool:
    """A drug or ligand treatment changes the chromatin; differentiation protocols do not count."""
    s = summary.lower()
    return "treated with" in s and "originated from" not in s


def _file_fields() -> str:
    fields = (
        "accession",
        "output_type",
        "file_format",
        "file_type",
        "file_size",
        "biological_replicates",
        "preferred_default",
        "cloud_metadata.url",
        "href",
    )
    return "".join(f"&field={f}" for f in fields)


def _file_row(f: dict) -> dict:
    url = (f.get("cloud_metadata") or {}).get("url") or (ENCODE + f.get("href", ""))
    return {
        "accession": f["accession"],
        "output_type": f.get("output_type"),
        "size": f.get("file_size"),
        "replicates": f.get("biological_replicates") or [],
        "url": url,
    }


def choose_histone_files(files: list[dict], analysis_files: set[str]) -> dict[str, dict | None]:
    """From an experiment's released GRCh38 files, the peak file and the signal file of its
    default analysis. Peaks: narrowPeak bed, replicated before pseudoreplicated, the
    portal's preferred default first. Signal: the fold-change bigWig pooling the most
    biological replicates, then the largest."""
    mine = [f for f in files if not analysis_files or f["accession"] in analysis_files]
    peaks = [
        f
        for f in mine
        if f.get("file_format") == "bed"
        and f.get("file_type") == "bed narrowPeak"
        and f.get("output_type") in PEAK_OUTPUTS
    ]
    peaks.sort(
        key=lambda f: (
            PEAK_OUTPUTS.index(f["output_type"]),
            not f.get("preferred_default"),
            -len(f.get("biological_replicates") or []),
            -(f.get("file_size") or 0),
        )
    )
    signal = [f for f in mine if f.get("file_format") == "bigWig" and f.get("output_type") == SIGNAL_OUTPUT]
    signal.sort(key=lambda f: (-len(f.get("biological_replicates") or []), -(f.get("file_size") or 0)))
    return {
        "peaks": _file_row(peaks[0]) if peaks else None,
        "signal": _file_row(signal[0]) if signal else None,
    }


def mark_options(cell_type: str, mark: str) -> list[dict]:
    """Every released GRCh38 Histone ChIP-seq experiment for one biosample and mark, with the
    peak and signal files of its default analysis (either may be None)."""
    q = (
        "/search/?type=Experiment&assay_title=Histone+ChIP-seq&status=released&limit=all&format=json"
        f"&biosample_ontology.term_name={urllib.parse.quote(cell_type)}"
        f"&target.label={urllib.parse.quote(mark)}"
        "&field=accession&field=assembly&field=biosample_summary&field=default_analysis"
        "&field=award.project&field=date_released"
    )
    exps = [
        e
        for e in _portal(q)["@graph"]
        if "GRCh38" in (e.get("assembly") or []) and not _treated(e.get("biosample_summary", ""))
    ]
    options = []
    for e in exps:
        analysis = e.get("default_analysis")
        afiles: set[str] = set()
        if analysis:
            a = _portal(f"{analysis}?format=json&frame=object")
            afiles = {p.strip("/").split("/")[-1] for p in a.get("files", [])}
        fq = (
            f"/search/?type=File&dataset=/experiments/{e['accession']}/&status=released&assembly=GRCh38"
            f"&limit=all&format=json{_file_fields()}"
        )
        chosen = choose_histone_files(_portal(fq)["@graph"], afiles)
        options.append(
            {
                "cell_type": cell_type,
                "mark": mark,
                "experiment": e["accession"],
                "biosample": e.get("biosample_summary", ""),
                "project": (e.get("award") or {}).get("project"),
                "released": e.get("date_released"),
                "analysis": analysis,
                "peaks": chosen["peaks"],
                "signal": chosen["signal"],
            }
        )
    return options


def pick_marks(options: dict[str, list[dict]]) -> dict[str, dict]:
    """One experiment per mark, from one biosample where possible.

    Marks from different donors or derivations of the same term (cardiac muscle from H7 or
    from RUES2, monocytes from two donors) do not describe one cell, so the biosample
    summary covering the most marks with complete files wins for every mark it covers;
    ties go to the larger total peak files. Inside a biosample: complete (peaks and signal)
    first, then the largest peak file."""

    def complete(o: dict) -> bool:
        return o["peaks"] is not None and o["signal"] is not None

    cover: dict[str, set[str]] = {}
    size: dict[str, int] = {}
    for mark, opts in options.items():
        for o in opts:
            if complete(o):
                cover.setdefault(o["biosample"], set()).add(mark)
                size[o["biosample"]] = size.get(o["biosample"], 0) + ((o["peaks"] or {}).get("size") or 0)
    ranked = sorted(cover, key=lambda b: (-len(cover[b]), -size[b], b))
    out: dict[str, dict] = {}
    for mark, opts in options.items():
        usable = [o for o in opts if o["peaks"] is not None or o["signal"] is not None]
        if not opts:
            out[mark] = {"mark": mark, "value": UNKNOWN, "reason": "no released GRCh38 experiment"}
            continue
        if not usable:
            out[mark] = {
                "mark": mark,
                "value": UNKNOWN,
                "reason": "experiments have no GRCh38 peaks or fold-change signal",
                "experiments": [o["experiment"] for o in opts],
            }
            continue

        def key(o: dict) -> tuple:
            rank = ranked.index(o["biosample"]) if o["biosample"] in ranked else len(ranked)
            return (not complete(o), rank, -((o["peaks"] or {}).get("size") or 0), o["experiment"])

        best = min(usable, key=key)
        out[mark] = {**best, "candidates": len(opts)}
    return out


def resolve_methylation(cell_type: str) -> dict:
    """The released GRCh38 WGBS experiment for one biosample (largest CpG bigBed), or why not."""
    q = (
        "/search/?type=Experiment&assay_title=WGBS&status=released&limit=all&format=json"
        f"&biosample_ontology.term_name={urllib.parse.quote(cell_type)}"
        "&field=accession&field=assembly&field=biosample_summary&field=award.project&field=date_released"
    )
    exps = [e for e in _portal(q)["@graph"] if "GRCh38" in (e.get("assembly") or [])]
    if not exps:
        rq = (
            "/search/?type=Experiment&assay_title=RRBS&status=released&limit=all&format=json"
            f"&biosample_ontology.term_name={urllib.parse.quote(cell_type)}&field=accession&field=assembly"
        )
        rrbs = _portal(rq)["@graph"]
        grch38 = [r for r in rrbs if "GRCh38" in (r.get("assembly") or [])]
        reason = "no released GRCh38 WGBS"
        reason += (
            f"; RRBS exists only on {sorted({a for r in rrbs for a in r.get('assembly') or []})}"
            if rrbs and not grch38
            else ("; no RRBS either" if not rrbs else "")
        )
        return {
            "cell_type": cell_type,
            "value": UNKNOWN,
            "reason": reason,
            "rrbs": [r["accession"] for r in rrbs],
        }
    best = None
    for e in exps:
        fq = (
            f"/search/?type=File&dataset=/experiments/{e['accession']}/&status=released&assembly=GRCh38"
            f"&file_format=bigBed&limit=all&format=json{_file_fields()}"
        )
        files = [f for f in _portal(fq)["@graph"] if f.get("output_type") == WGBS_OUTPUT]
        files.sort(key=lambda f: (not f.get("preferred_default"), -(f.get("file_size") or 0)))
        if files and (best is None or (files[0].get("file_size") or 0) > (best[1].get("file_size") or 0)):
            best = (e, files[0])
    if best is None:
        return {"cell_type": cell_type, "value": UNKNOWN, "reason": "WGBS without a GRCh38 CpG bigBed"}
    e, f = best
    return {
        "cell_type": cell_type,
        "experiment": e["accession"],
        "biosample": e.get("biosample_summary", ""),
        "project": (e.get("award") or {}).get("project"),
        "released": e.get("date_released"),
        "cpg": _file_row(f),
        "candidates": len(exps),
    }


def build_manifest(cell_types: Iterable[str] = CELL_TYPES, progress=None) -> dict:
    """Resolve every (biosample, mark) and every biosample's WGBS; saved as the pinned manifest."""
    from genomeos.results import save_result

    cells: dict[str, dict] = {}
    for cell in cell_types:
        row: dict[str, Any] = {"dnase": dnase_accession(cell)}
        options = {}
        for mark in MARKS:
            options[mark] = mark_options(cell, mark)
            if progress:
                progress(cell, mark)
        row["marks"] = {m: {"cell_type": cell, **e} for m, e in pick_marks(options).items()}
        row["methylation"] = resolve_methylation(cell)
        cells[cell] = row
    manifest = {
        "source": ENCODE,
        "rule": {
            "experiment": "released, GRCh38, no drug treatment; per cell type the biosample covering the "
            "most marks, then complete files (peaks and signal), then the largest peak file",
            "peaks": "default analysis, replicated narrowPeak before pseudoreplicated",
            "signal": "default analysis, fold change over control, most pooled replicates",
            "methylation": "WGBS 'methylation state at CpG' bigBed, the largest over the biosample's "
            "experiments",
        },
        "cell_types": cells,
        "coverage": coverage_table(cells),
    }
    save_result(MANIFEST, manifest)
    return manifest


def load_manifest() -> dict | None:
    from genomeos.results import load_result

    return load_result(MANIFEST)


def dnase_accession(cell_type: str) -> str | None:
    """The DNase file the reader streamed for this biosample, from its local peak cache header."""
    for p in sorted(RESULTS.glob(f"dnase_{slug(cell_type)}_chr*.bed.gz")):
        with gzip.open(p, "rt") as fh:
            head = fh.readline()
        if head.startswith("#"):
            return head[1:].split()[0]
    return None


def coverage_table(cells: dict[str, dict]) -> dict:
    """Which of the biosamples have which experiment: the honest size of the layer."""
    rows = {}
    for cell, row in cells.items():
        marks = {m: ("value" not in r) for m, r in row["marks"].items()}
        rows[cell] = {
            "dnase": bool(row.get("dnase")),
            **marks,
            "methylation": "value" not in row["methylation"],
        }
    n = len(rows)
    complete = [c for c, r in rows.items() if all(r.values())]
    return {
        "rows": rows,
        "biosamples": n,
        "with_all_five_marks": sum(1 for r in rows.values() if all(r[m] for m in MARKS)),
        "with_methylation": sum(1 for r in rows.values() if r["methylation"]),
        "complete": complete,
        "missing": {c: [k for k, v in r.items() if not v] for c, r in rows.items() if not all(r.values())},
    }


# ------------------------------------------------------------------------------------------
# Peaks: streamed once per (biosample, mark), rows kept per chromosome
# ------------------------------------------------------------------------------------------


def peaks_path(cell_type: str, mark: str, chrom: str) -> Path:
    return CACHE / "peaks" / f"{slug(cell_type)}_{mark}_{chrom}.bed.gz"


def fetch_mark_peaks(entry: dict, chroms: Iterable[str] = CHROMS, timeout: int = 300) -> dict[str, int]:
    """Stream one narrowPeak file and keep (start, end, signalValue) per chromosome."""
    cell, mark, info = entry["cell_type"], entry["mark"], entry["peaks"]
    wanted = set(chroms)
    kept: dict[str, list[tuple[int, int, float]]] = {c: [] for c in wanted}
    req = urllib.request.Request(info["url"], headers={"User-Agent": USER_AGENT})
    with (
        urllib.request.urlopen(req, timeout=timeout) as resp,  # noqa: S310
        gzip.open(io.BufferedReader(resp, 1 << 20), "rt") as fh,
    ):
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or f[0] not in kept:
                continue
            kept[f[0]].append((int(f[1]), int(f[2]), float(f[6])))
    for c, rows in kept.items():
        rows.sort()
        p = peaks_path(cell, mark, c)
        p.parent.mkdir(parents=True, exist_ok=True)
        part = p.with_name(p.name + ".part")
        with gzip.open(part, "wt") as out:
            out.write(
                f"# {info['accession']} ({entry['experiment']}) {cell} {mark} {info['output_type']}, {c}\n"
            )
            for s, e, v in rows:
                out.write(f"{s}\t{e}\t{v:.3f}\n")
        part.replace(p)  # atomic: an interrupted run never leaves a half file that reads as cached
    return {c: len(v) for c, v in kept.items()}


def load_mark_peaks(cell_type: str, mark: str, chrom: str) -> list[tuple[int, int, float]] | None:
    """Peaks of one mark on one chromosome; None when never fetched (not the same as no peaks)."""
    p = peaks_path(cell_type, mark, chrom)
    if not p.exists():
        return None
    out = []
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            s, e, v = line.rstrip("\n").split("\t")
            out.append((int(s), int(e), float(v)))
    return out


# ------------------------------------------------------------------------------------------
# Signal: fold change over control, 200 bp bin means per chromosome
# ------------------------------------------------------------------------------------------


def signal_path(cell_type: str, mark: str, chrom: str) -> Path:
    return CACHE / "signal" / f"{slug(cell_type)}_{mark}_{chrom}.f32.gz"


def fetch_signal_profile(entry: dict, chrom: str) -> dict:
    """Read one chromosome of a fold-change bigWig by range into BIN-bp means (NaN: no value)."""
    bw = _bw.BigWig(entry["signal"]["url"])
    t0 = time.time()
    try:
        if chrom not in bw.chroms:
            prof = array("f")
        else:
            length = bw.chroms[chrom][1]
            bins = [(s, min(s + BIN, length)) for s in range(0, length, BIN)]
            stats = bw.summarise(chrom, bins, float("inf"))
            prof = array("f", (st.mean if st.bases else math.nan for st in stats))
        cost = {"requests": bw.src.requests, "mb": round(bw.src.bytes_fetched / 1e6, 2)}
    finally:
        bw.close()
    p = signal_path(entry["cell_type"], entry["mark"], chrom)
    p.parent.mkdir(parents=True, exist_ok=True)
    part = p.with_name(p.name + ".part")
    with gzip.open(part, "wb") as fh:
        fh.write(prof.tobytes())
    part.replace(p)
    return {**cost, "seconds": round(time.time() - t0, 1), "bins": len(prof)}


def signal_coverage(cache: Path | None = None) -> dict[str, dict[str, list[str]]]:
    """Which chromosomes carry a fold-change profile, per cell type and mark (the rest: peaks only)."""
    root = (cache or CACHE) / "signal"
    out: dict[str, dict[str, list[str]]] = {}
    for cell in CELL_TYPES:
        for mark in MARKS:
            prefix = f"{slug(cell)}_{mark}_"
            chroms = sorted(
                (p.name[len(prefix) : -len(".f32.gz")] for p in root.glob(prefix + "chr*.f32.gz")),
                key=lambda c: CHROMS.index(c) if c in CHROMS else 99,
            )
            out.setdefault(cell, {})[mark] = chroms
    return out


def load_signal_profile(cell_type: str, mark: str, chrom: str) -> array | None:
    p = signal_path(cell_type, mark, chrom)
    if not p.exists():
        return None
    prof = array("f")
    with gzip.open(p, "rb") as fh:
        prof.frombytes(fh.read())
    return prof


def profile_mean(prof: array, start: int, end: int) -> float | None:
    """Mean of the bin means a [start, end) interval touches; None when none carries a value."""
    lo, hi = start // BIN, (max(end, start + 1) - 1) // BIN + 1
    vals = [v for v in prof[lo:hi] if v == v]
    return sum(vals) / len(vals) if vals else None


# ------------------------------------------------------------------------------------------
# Methylation: a bigBed range reader on the bigWig machinery, CpG calls into bins
# ------------------------------------------------------------------------------------------


class BigBed(_bw.BigWig):
    """A bigBed opened by range: the header, chromosome tree and R-tree are the bigWig ones;
    only the magic and the data sections (bed rows instead of values) differ."""

    def __init__(self, source: str | Path):  # noqa: D107 - mirrors BigWig.__init__
        self.src = _bw._Source(source)
        h = self.src.read(0, 64)
        magic, version, zoom, ct, fd, fi, fc, _dfc, _asql, _ts, ubs, _res = struct.unpack("<IHHQQQHHQQIQ", h)
        if magic != BIGBED_MAGIC:
            raise ValueError(f"not a little-endian bigBed: magic {magic:#x}")
        self.header = _bw.Header(version, zoom, ct, fd, fi, ubs)
        self.field_count = fc
        self.chroms = {}
        self._read_chrom_tree()
        (_m, self._rtree_block_size, _n, _sc, _sb, _ec, _eb, _eo, self._items_per_slot, _r) = struct.unpack(
            "<IIQIIIIQII", self.src.read(fi, 48)
        )
        self._node_cache = {}

    def rows(self, chrom: str, intervals: list[tuple[int, int]]) -> Iterator[tuple[int, int, str]]:
        """Every (start, end, rest-of-line) of `chrom` in the data blocks overlapping the intervals."""
        if chrom not in self.chroms:
            return
        cid = self.chroms[chrom][0]
        for group in _bw.coalesce(self.leaf_items(chrom, intervals)):
            base = group[0].offset
            blob = self.src.read(base, group[-1].offset + group[-1].size - base)
            for it in group:
                raw = blob[it.offset - base : it.offset - base + it.size]
                if self.header.uncompress_buf_size:
                    raw = zlib.decompress(raw)
                yield from parse_bed_block(raw, cid)


def parse_bed_block(raw: bytes, chrom_id: int) -> Iterator[tuple[int, int, str]]:
    """A decompressed bigBed data block: (chromId, start, end, rest\\0) records."""
    p, n = 0, len(raw)
    while p + 12 <= n:
        cid, s, e = struct.unpack_from("<III", raw, p)
        z = raw.index(b"\0", p + 12)
        if cid == chrom_id:
            yield s, e, raw[p + 12 : z].decode("ascii", "replace")
        p = z + 1


def bedmethyl_call(rest: str) -> tuple[int, float] | None:
    """(coverage, fraction) from the columns after chrom/start/end of an ENCODE bedMethyl row:
    name, score, strand, thickStart, thickEnd, colour, coverage, percent, ..."""
    f = rest.split("\t")
    if len(f) < 8:
        return None
    try:
        return int(f[6]), float(f[7]) / 100.0
    except ValueError:
        return None


# per bin: CpG calls, calls with >= MIN_COVERAGE reads, sum of their fractions x 1000, their reads
METH_COLUMNS = ("calls", "covered", "fraction_sum_permille", "reads")


def methylation_path(cell_type: str, chrom: str) -> Path:
    return CACHE / "methylation" / f"{slug(cell_type)}_{chrom}.u32.gz"


def bin_calls(calls: Iterable[tuple[int, int, float]], length: int) -> list[array]:
    """(position, coverage, fraction) calls into BIN-bp columns (METH_COLUMNS)."""
    nb = length // BIN + 1
    cols = [array("I", bytes(4 * nb)) for _ in METH_COLUMNS]
    calls_c, cov_c, frac_c, reads_c = cols
    for pos, cov, frac in calls:
        b = pos // BIN
        if b >= nb:
            continue
        calls_c[b] += 1
        if cov >= MIN_COVERAGE:
            cov_c[b] += 1
            frac_c[b] += round(frac * 1000)
            reads_c[b] += cov
    return cols


def fetch_methylation_profile(entry: dict, chrom: str) -> dict:
    """Read one chromosome of the CpG bigBed by range and keep the binned columns."""
    bb = BigBed(entry["cpg"]["url"])
    t0 = time.time()
    try:
        if chrom not in bb.chroms:
            cols = [array("I") for _ in METH_COLUMNS]
            n = 0
        else:
            length = bb.chroms[chrom][1]
            n = 0

            def calls() -> Iterator[tuple[int, int, float]]:
                nonlocal n
                for s, _e, rest in bb.rows(chrom, [(0, length)]):
                    c = bedmethyl_call(rest)
                    if c:
                        n += 1
                        yield s, c[0], c[1]

            cols = bin_calls(calls(), length)
        cost = {"requests": bb.src.requests, "mb": round(bb.src.bytes_fetched / 1e6, 2)}
    finally:
        bb.close()
    p = methylation_path(entry["cell_type"], chrom)
    p.parent.mkdir(parents=True, exist_ok=True)
    part = p.with_name(p.name + ".part")
    with gzip.open(part, "wb") as fh:
        fh.write(struct.pack("<I", len(cols[0])))
        for c in cols:
            fh.write(c.tobytes())
    part.replace(p)
    return {**cost, "seconds": round(time.time() - t0, 1), "calls": n}


def load_methylation_profile(cell_type: str, chrom: str) -> list[array] | None:
    p = methylation_path(cell_type, chrom)
    if not p.exists():
        return None
    with gzip.open(p, "rb") as fh:
        raw = fh.read()
    (nb,) = struct.unpack_from("<I", raw, 0)
    cols = []
    for i in range(len(METH_COLUMNS)):
        a = array("I")
        a.frombytes(raw[4 + i * 4 * nb : 4 + (i + 1) * 4 * nb])
        cols.append(a)
    return cols


def methylation_over(cols: list[array], start: int, end: int) -> dict:
    """Methylation of the bins [start, end) touches: mean per-call fraction over calls with at
    least MIN_COVERAGE reads, how many calls there are and their mean coverage."""
    lo, hi = start // BIN, (max(end, start + 1) - 1) // BIN + 1
    calls = sum(cols[0][lo:hi])
    covered = sum(cols[1][lo:hi])
    frac = sum(cols[2][lo:hi])
    reads = sum(cols[3][lo:hi])
    return {
        "fraction": round(frac / covered / 1000, 4) if covered else None,
        "cpg_calls": calls,
        "cpg_calls_covered": covered,
        "mean_coverage": round(reads / covered, 1) if covered else None,
    }


def cpg_sites(seq: str) -> int:
    """CpG dinucleotides in a sequence (each is two strand calls in a bedMethyl file)."""
    return seq.upper().count("CG")


# ------------------------------------------------------------------------------------------
# The record: one locus, one cell type, every field with evidence and confidence
# ------------------------------------------------------------------------------------------


class Layer:
    """Everything cached for a set of cell types on one chromosome, loaded once."""

    def __init__(
        self,
        chrom: str,
        cell_types: Iterable[str] = CELL_TYPES,
        manifest: dict | None = None,
        signal: bool = True,
        methylation: bool = True,
    ):
        self.chrom = chrom
        self.manifest = manifest if manifest is not None else (load_manifest() or {"cell_types": {}})
        self.cells = list(cell_types)
        self.dnase: dict[str, reader.PeakIndex | None] = {}
        self.peaks: dict[tuple[str, str], reader.PeakIndex | None] = {}
        self.signal: dict[tuple[str, str], array | None] = {}
        self.meth: dict[str, list[array] | None] = {}
        for cell in self.cells:
            d = reader.load_peaks(cell, chrom)
            self.dnase[cell] = reader.PeakIndex(d) if d else None
            for mark in MARKS:
                pk = load_mark_peaks(cell, mark, chrom)
                self.peaks[(cell, mark)] = reader.PeakIndex(pk) if pk is not None else None
                self.signal[(cell, mark)] = load_signal_profile(cell, mark, chrom) if signal else None
            self.meth[cell] = load_methylation_profile(cell, chrom) if methylation else None

    def entry(self, cell: str, mark: str | None = None) -> dict:
        row = self.manifest.get("cell_types", {}).get(cell, {})
        if mark is None:
            return row.get("methylation") or {"value": UNKNOWN, "reason": "not in the manifest"}
        return (row.get("marks") or {}).get(mark) or {"value": UNKNOWN, "reason": "not in the manifest"}

    def record(self, start: int, end: int, cell: str, name: str = "", kind: str = "locus") -> dict:
        return {
            "locus": {"chrom": self.chrom, "start": start, "end": end, "name": name, "kind": kind},
            "cell_type": cell,
            "openness": self._openness(cell, start, end),
            "marks": {m: self._mark(cell, m, start, end) for m in MARKS},
            "methylation": self._methylation(cell, start, end),
            "state": chromatin_state_field(self._marks_called(cell, start, end)),
        }

    def _openness(self, cell: str, start: int, end: int) -> dict:
        idx = self.dnase.get(cell)
        if idx is None:
            return {
                "value": UNKNOWN,
                "reason": f"no DNase peaks read for {cell} on {self.chrom} (genomeos reader)",
            }
        hits = idx.overlapping(start, end)
        return {
            "value": "open" if hits else "closed",
            "dnase_signal": round(max((v for _, _, v in hits), default=0.0), 2),
            "evidence": reader.EVIDENCE,
            "confidence": 0.8,
            "source": dnase_accession(cell),
        }

    def _marks_called(self, cell: str, start: int, end: int) -> dict[str, bool | None]:
        out: dict[str, bool | None] = {}
        for m in MARKS:
            idx = self.peaks.get((cell, m))
            out[m] = None if idx is None else bool(idx.overlapping(start, end))
        return out

    def _mark(self, cell: str, mark: str, start: int, end: int) -> dict:
        e = self.entry(cell, mark)
        if "value" in e:
            return {"value": UNKNOWN, "reason": e.get("reason", "no experiment")}
        idx = self.peaks.get((cell, mark))
        prof = self.signal.get((cell, mark))
        if idx is None and prof is None:
            return {"value": UNKNOWN, "reason": "experiment exists; this chromosome not read yet"}
        out: dict[str, Any] = {
            "mark": MARKS[mark],
            "experiment": e["experiment"],
            "biosample": e.get("biosample"),
        }
        if idx is not None:
            hits = idx.overlapping(start, end)
            out["peak"] = bool(hits)
            out["peak_signal"] = round(max((v for _, _, v in hits), default=0.0), 3)
            out["peak_file"] = (e.get("peaks") or {}).get("accession")
        else:
            out["peak"] = UNKNOWN
        if prof is not None:
            fc = profile_mean(prof, start, end)
            out["fold_change"] = round(fc, 3) if fc is not None else UNKNOWN
            out["signal_file"] = (e.get("signal") or {}).get("accession")
        else:
            out["fold_change"] = UNKNOWN
            out["fold_change_reason"] = "fold-change profile not read for this chromosome; peaks only"
        replicated = (e.get("peaks") or {}).get("output_type") == "replicated peaks"
        pooled = len((e.get("signal") or {}).get("replicates") or []) >= 2
        out["evidence"] = EVIDENCE_MARK
        out["confidence"] = round(0.9 if replicated and pooled else 0.8 if replicated or pooled else 0.7, 2)
        return out

    def _methylation(self, cell: str, start: int, end: int) -> dict:
        e = self.entry(cell)
        if "value" in e:
            return {"value": UNKNOWN, "reason": e.get("reason", "no experiment")}
        cols = self.meth.get(cell)
        if cols is None:
            return {"value": UNKNOWN, "reason": "WGBS exists; this chromosome not read yet"}
        m = methylation_over(cols, start, end)
        if not m["cpg_calls_covered"]:
            return {
                "value": UNKNOWN,
                "reason": f"no CpG call with {MIN_COVERAGE}+ reads here",
                **m,
                "experiment": e["experiment"],
            }
        conf = 0.9 if m["cpg_calls_covered"] >= 10 and (m["mean_coverage"] or 0) >= 10 else 0.6
        return {
            **m,
            "experiment": e["experiment"],
            "biosample": e.get("biosample"),
            "file": e["cpg"]["accession"],
            "evidence": EVIDENCE_METHYLATION + f"; calls with {MIN_COVERAGE}+ reads, {BIN} bp resolution",
            "confidence": conf,
        }


def chromatin_state(called: dict[str, bool | None]) -> str:
    """A state from the five marks' peak calls at a locus; UNKNOWN when a mark the rule needs is
    unmeasured. The rules are the textbook combinations, applied in order."""
    k4me3, k27ac, k4me1 = called.get("H3K4me3"), called.get("H3K27ac"), called.get("H3K4me1")
    k27me3, k9me3 = called.get("H3K27me3"), called.get("H3K9me3")
    if k4me3 and k27me3:
        return "bivalent"
    if k4me3:
        return "active promoter" if k27ac else "promoter mark, not acetylated"
    if k27ac:
        return "active enhancer" if k4me1 else "acetylated element"
    if k4me1:
        return "poised enhancer" if k27me3 else "primed enhancer"
    if k27me3:
        return "Polycomb-repressed"
    if k9me3:
        return "heterochromatin"
    if any(v is None for v in called.values()):
        return UNKNOWN
    return "no mark"


STATES = (
    "bivalent",
    "active promoter",
    "promoter mark, not acetylated",
    "active enhancer",
    "acetylated element",
    "poised enhancer",
    "primed enhancer",
    "Polycomb-repressed",
    "heterochromatin",
    "no mark",
    UNKNOWN,
)


def chromatin_state_field(called: dict[str, bool | None]) -> dict:
    s = chromatin_state(called)
    return {
        "value": s,
        "marks_called": [m for m, v in called.items() if v],
        "marks_unmeasured": [m for m, v in called.items() if v is None],
        "evidence": EVIDENCE_STATE,
        "confidence": 0.5 if s != UNKNOWN else 0.0,
    }


def epigenome_at(
    chrom: str,
    start: int,
    end: int,
    cell_types: Iterable[str] = CELL_TYPES,
    name: str = "",
    kind: str = "locus",
    layer: Layer | None = None,
) -> list[dict]:
    """The epigenome records of one locus: one per cell type, from what is cached locally.

    Fields a biosample has no experiment for, or that were never read for this chromosome,
    come back as UNKNOWN with the reason; nothing is fetched here (`ensure` does that)."""
    cells = list(cell_types)
    layer = layer or Layer(chrom, cells)
    return [layer.record(start, end, cell, name=name, kind=kind) for cell in cells]


def ensure(
    chrom: str,
    cell_types: Iterable[str] = CELL_TYPES,
    signal: bool = True,
    methylation: bool = True,
    log=print,
) -> dict:
    """Fetch what is missing for one chromosome: peaks (all chromosomes at once), signal
    profiles and methylation bins. Returns what it cost."""
    manifest = load_manifest() or build_manifest(cell_types)
    cost: dict[str, Any] = {"peaks": 0, "signal_mb": 0.0, "methylation_mb": 0.0, "seconds": 0.0}
    t0 = time.time()
    for cell in cell_types:
        row = manifest["cell_types"].get(cell) or {}
        for mark, e in (row.get("marks") or {}).items():
            if "value" in e:
                continue
            if e.get("peaks") and load_mark_peaks(cell, mark, chrom) is None:
                kept = fetch_mark_peaks(e)
                cost["peaks"] += 1
                log(f"  {cell} {mark}: {e['peaks']['accession']} peaks, {kept.get(chrom, 0):,} on {chrom}")
            if signal and e.get("signal") and load_signal_profile(cell, mark, chrom) is None:
                c = fetch_signal_profile(e, chrom)
                cost["signal_mb"] += c["mb"]
                acc = e["signal"]["accession"]
                log(f"  {cell} {mark}: {acc} fold change on {chrom}, {c['mb']} MB {c['seconds']} s")
        m = row.get("methylation") or {}
        if methylation and "cpg" in m and load_methylation_profile(cell, chrom) is None:
            c = fetch_methylation_profile(m, chrom)
            cost["methylation_mb"] += c["mb"]
            acc = m["cpg"]["accession"]
            log(f"  {cell} WGBS: {acc} on {chrom}, {c['calls']:,} CpG calls, {c['mb']} MB {c['seconds']} s")
    cost["seconds"] = round(time.time() - t0, 1)
    cost["signal_mb"] = round(cost["signal_mb"], 1)
    cost["methylation_mb"] = round(cost["methylation_mb"], 1)
    return cost


# ------------------------------------------------------------------------------------------
# The per-chromosome summary
# ------------------------------------------------------------------------------------------


def _tss(g) -> int:
    from genomeos.coords import Strand

    return g.locus.end - 1 if g.locus.strand is Strand.MINUS else g.locus.start


def _share(n: int, d: int) -> float | None:
    return round(n / d, 4) if d else None


def _median(xs: list[float]) -> float | None:
    xs = sorted(xs)
    if not xs:
        return None
    m = len(xs) // 2
    return round(xs[m] if len(xs) % 2 else (xs[m - 1] + xs[m]) / 2, 4)


def summarise_chromosome(
    chrom: str,
    annotation,
    ccres: list,
    budget: dict | None = None,
    cell_types: Iterable[str] = CELL_TYPES,
    manifest: dict | None = None,
) -> dict:
    """What the layer says about one chromosome, per cell type: marks on promoters (read and
    silent by the DNase reader), marks on the registry's elements by class, promoter
    methylation, and methylation and coverage per budget tier of the UNKNOWN space."""
    cells = list(cell_types)
    layer = Layer(chrom, cells, manifest=manifest, signal=False)
    genes = sorted(
        (g for g in annotation.genes.values() if g.locus.chrom == chrom and g.type == "protein_coding"),
        key=_tss,
    )
    out_cells: dict[str, Any] = {}
    for cell in cells:
        dn = layer.dnase.get(cell)
        measured = [m for m in MARKS if layer.peaks.get((cell, m)) is not None]
        row: dict[str, Any] = {
            "marks_measured": measured,
            "marks_unknown": {
                m: layer._mark(cell, m, 0, 1).get("reason") for m in MARKS if m not in measured
            },
            "dnase": dn is not None,
            "methylation": layer.meth.get(cell) is not None,
        }
        row["peaks"] = {
            m: {
                "n": len(layer.peaks[(cell, m)].peaks),
                "bp": sum(e - s for s, e, _ in layer.peaks[(cell, m)].peaks),
            }
            for m in measured
        }
        # promoters, split by the reader's call
        prom: dict[str, dict[str, Any]] = {"read": {}, "silent": {}, "unread": {}}
        meth_by: dict[str, list[float]] = {"read": [], "silent": []}
        for g in genes:
            t = _tss(g)
            lo, hi = max(0, t - PROMOTER_WINDOW), t + PROMOTER_WINDOW
            group = "unread" if dn is None else ("read" if dn.overlapping(lo, hi) else "silent")
            called = layer._marks_called(cell, lo, hi)
            bucket = prom[group]
            bucket["n"] = bucket.get("n", 0) + 1
            for m, v in called.items():
                if v:
                    bucket[m] = bucket.get(m, 0) + 1
            st = chromatin_state(called)
            bucket.setdefault("states", {})[st] = bucket.setdefault("states", {}).get(st, 0) + 1
            cols = layer.meth.get(cell)
            if cols is not None and group in meth_by:
                mm = methylation_over(cols, lo, hi)
                if mm["cpg_calls_covered"] >= 4:
                    meth_by[group].append(mm["fraction"])
        prom = {k: v for k, v in prom.items() if v.get("n")}
        for bucket in prom.values():
            for m in measured:
                bucket[m] = {"n": bucket.get(m, 0), "share": _share(bucket.get(m, 0), bucket["n"])}
        row["promoters"] = prom
        if layer.meth.get(cell) is not None:
            row["promoter_methylation_median"] = {k: _median(v) for k, v in meth_by.items()}
            row["promoter_methylation_n"] = {k: len(v) for k, v in meth_by.items()}
        # registry elements by class
        by_cls: dict[str, dict[str, Any]] = {}
        for c in ccres:
            called = layer._marks_called(cell, c.start, c.end)
            b = by_cls.setdefault(c.cls, {"n": 0, "open": 0, "states": {}})
            b["n"] += 1
            if dn is not None and dn.overlapping(c.start, c.end):
                b["open"] += 1
            for m, v in called.items():
                if v:
                    b[m] = b.get(m, 0) + 1
            st = chromatin_state(called)
            b["states"][st] = b["states"].get(st, 0) + 1
        for b in by_cls.values():
            for m in measured:
                b[m] = {"n": b.get(m, 0), "share": _share(b.get(m, 0), b["n"])}
            b["open_share"] = _share(b["open"], b["n"]) if dn is not None else None
        row["registry_classes"] = by_cls
        # methylation per budget tier
        cols = layer.meth.get(cell)
        if cols is not None:
            total_calls = sum(cols[0])
            total_cov = sum(cols[1])
            row["cpg_calls"] = total_calls
            row["cpg_calls_covered"] = total_cov
            if budget:
                tiers: dict[str, dict[str, Any]] = {}
                for blk in budget.get("blocks", []):
                    tier = blk["guess"]["tier"]
                    m = methylation_over(cols, blk["start"], blk["end"])
                    t = tiers.setdefault(tier, {"blocks": 0, "bp": 0, "calls": 0, "covered": 0, "_fsum": 0.0})
                    t["blocks"] += 1
                    t["bp"] += blk["end"] - blk["start"]
                    t["calls"] += m["cpg_calls"]
                    t["covered"] += m["cpg_calls_covered"]
                    if m["fraction"] is not None:
                        t["_fsum"] += m["fraction"] * m["cpg_calls_covered"]
                for t in tiers.values():
                    t["fraction"] = round(t.pop("_fsum") / t["covered"], 4) if t["covered"] else None
                row["methylation_by_tier"] = tiers
        out_cells[cell] = row
    return {
        "chrom": chrom,
        "coding_genes": len(genes),
        "registry_elements": len(ccres),
        "cell_types": out_cells,
        "evidence": {
            "marks": EVIDENCE_MARK + ", replicated peaks overlapping the locus",
            "promoters": "TSS +/- 1 kb; read or silent by the DNase reader (inferred)",
            "methylation": EVIDENCE_METHYLATION + f", calls with {MIN_COVERAGE}+ reads, promoters with 4+",
            "state": EVIDENCE_STATE,
        },
    }


# ------------------------------------------------------------------------------------------
# Small statistics, standard library only (used by scripts/epigenome.py)
# ------------------------------------------------------------------------------------------


def ranks(xs: list[float]) -> list[float]:
    """Average ranks (1-based), ties shared."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def auc(scores: list[float], labels: list[int]) -> float | None:
    """Area under the ROC curve (Mann-Whitney), ties counted half."""
    pos = sum(labels)
    neg = len(labels) - pos
    if not pos or not neg:
        return None
    r = ranks(scores)
    return (sum(ri for ri, y in zip(r, labels, strict=True) if y) - pos * (pos + 1) / 2) / (pos * neg)


def spearman(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3:
        return None
    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else None


def solve(a: list[list[float]], b: list[float]) -> list[float]:
    """Gaussian elimination with partial pivoting (small dense systems)."""
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[p] = m[p], m[c]
        if abs(m[c][c]) < 1e-12:
            continue
        for r in range(n):
            if r != c and m[r][c]:
                f = m[r][c] / m[c][c]
                for k in range(c, n + 1):
                    m[r][k] -= f * m[c][k]
    return [m[i][n] / m[i][i] if abs(m[i][i]) >= 1e-12 else 0.0 for i in range(n)]


def ridge_fit(x: list[list[float]], y: list[float], lam: float = 1.0) -> list[float]:
    """Least squares with a small ridge on every coefficient but the intercept (column 0)."""
    k = len(x[0])
    xtx = [[0.0] * k for _ in range(k)]
    xty = [0.0] * k
    for row, yi in zip(x, y, strict=True):
        for i in range(k):
            ri = row[i]
            if ri == 0.0:
                continue
            xty[i] += ri * yi
            xi = xtx[i]
            for j in range(i, k):
                xi[j] += ri * row[j]
    for i in range(k):
        for j in range(i):
            xtx[i][j] = xtx[j][i]
        if i:
            xtx[i][i] += lam
    return solve(xtx, xty)


def predict(w: list[float], x: list[list[float]]) -> list[float]:
    return [sum(a * b for a, b in zip(w, row, strict=True)) for row in x]


def cross_validated(
    x: list[list[float]], y: list[float], groups: list[int], folds: int = 5, lam: float = 1.0
) -> list[float]:
    """Out-of-fold predictions of a ridge fit; rows sharing a group stay in one fold."""
    pred = [0.0] * len(y)
    for f in range(folds):
        tr = [i for i, g in enumerate(groups) if g % folds != f]
        te = [i for i, g in enumerate(groups) if g % folds == f]
        if not tr or not te:
            continue
        w = ridge_fit([x[i] for i in tr], [y[i] for i in tr], lam)
        for i, p in zip(te, predict(w, [x[i] for i in te]), strict=True):
            pred[i] = p
    return pred


def standardise(x: list[list[float]]) -> list[list[float]]:
    """Columns to zero mean and unit variance, with an intercept column prepended."""
    k = len(x[0]) if x else 0
    n = len(x)
    means = [sum(r[j] for r in x) / n for j in range(k)]
    sds = [math.sqrt(sum((r[j] - means[j]) ** 2 for r in x) / n) or 1.0 for j in range(k)]
    return [[1.0, *((r[j] - means[j]) / sds[j] for j in range(k))] for r in x]


def gc_cpg(seq: str) -> tuple[float, float]:
    """GC fraction over called bases and CpG dinucleotides per 100 called bases."""
    s = seq.upper()
    acgt = sum(s.count(b) for b in "ACGT")
    if not acgt:
        return math.nan, math.nan
    return (s.count("G") + s.count("C")) / acgt, 100 * s.count("CG") / acgt


def stratum(
    gc: float, cpg: float, gc_step: float = 0.05, cpg_edges: tuple[float, ...] = (0.5, 1, 2, 4, 8)
) -> str:
    """A GC by CpG-density stratum label."""
    return f"gc{int(gc / gc_step)}_cpg{bisect.bisect_right(list(cpg_edges), cpg)}"


__all__ = [
    "CELL_TYPES",
    "MARKS",
    "UNKNOWN",
    "BigBed",
    "Layer",
    "build_manifest",
    "chromatin_state",
    "ensure",
    "epigenome_at",
    "load_manifest",
    "summarise_chromosome",
]
