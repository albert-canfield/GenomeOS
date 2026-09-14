# SPDX-License-Identifier: AGPL-3.0-or-later
"""Many human genomes at once: what is fixed, what holds a value, what cannot be placed.

The project's human axis so far was gnomAD's Gnocchi (attribution/variation.py), one Z score per
kilobase from 76,156 genomes, with SNVs only and kilobases left unscored. This module reads the
genomes themselves. The Human Pangenome Reference Consortium's first release aligned to hg38 by
Cactus (UCSC `hprc90way`: hg38, T2T-CHM13 and 88 haplotype assemblies of 44 people, Liao et al.
2023) says, for every hg38 base, in how many of the 89 other assemblies it is present, deleted,
replaced or missing, and which allele each one carries. The alignment is 250 GB of MAF text at
UCSC (3.1 GB for chromosome 21) and never enters the repository: UCSC's table of block offsets is
read through the REST API, the byte ranges are streamed, and each block is distilled on the spot
into three small things that are kept under `data/knowledge/human_panel/<chrom>` (git-ignored):
the block's span and the status of every assembly in it, the columns where some assembly differs
from hg38 with the allele of each assembly, and the bases each assembly inserts. Everything below
is computed from that store.

Two readings come out of it.

*Blocks* (the budget's UNKNOWN blocks, with GENCODE coding exons and the neutral tier as the
controls) get the numbers the question asks for, presence, identity, the share of bases touched by
a deletion, replacement or insertion, and the density of common variable sites against a
background matched for GC, and a class: core, variable, polymorphic, lineage-restricted, or
unplaced when the panel does not see the sequence.

*Units* (fixed windows tiled over the same space) get a value domain: for each assembly the
combination of states it carries at the unit's recurring variable sites and inserted lengths, how
many distinct values recur, the commonest value's share, the entropy and the effective number of
values. A unit is fixed when one value holds 95% of the assemblies, storage when a handful of
recurring values covers them (Albert's multiple-select column, 2026-09-13), hypervariable when
most assemblies carry a value of their own, and unplaced when too few assemblies inform it.
Frequencies beyond 89 haplotypes come from population sources read by range as bigBed: gnomAD
v4.1.1 genomes for substitutions and short indels, gnomAD v4.1 structural variants, and the
TRExplorer catalogue for tandem-repeat allele sizes; each carries its own sample size and is never
mixed with the panel's counts. Design and numbers: docs/ATTRIBUTION.md, "Many human genomes".
"""

from __future__ import annotations

import contextlib
import functools
import gzip
import http.client
import json
import math
import re
import struct
import threading
import time
import urllib.parse
import zlib
from bisect import bisect_left, bisect_right
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.attribution.bigwig import BigWig, Header, _Source, coalesce
from genomeos.attribution.constraint import _api

# -- sources -------------------------------------------------------------------------------------
MAF_TRACK = "hprc90way"
MAF_URL = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/hprc/cactus90way/{chrom}.maf"
GBDB = "https://hgdownload.soe.ucsc.edu/gbdb/hg38"
REFERENCE = "hg38"
BIGBEDS = {
    "hprc_dup": f"{GBDB}/hprcArrV1/hprcArrDupV1.bb",
    "hprc_inv": f"{GBDB}/hprcArrV1/hprcArrInvV1.bb",
    "hprc_del": f"{GBDB}/hprcArrV1/hprcDeletionsV1.bb",
    "hprc_ins": f"{GBDB}/hprcArrV1/hprcInsertsV1.bb",
    "hprc_double": f"{GBDB}/hprcArrV1/hprcDoubleV1.bb",
    "hprc_v21_sv": f"{GBDB}/lrSv/hprc2v21.bb",
    "gnomad_snv": f"{GBDB}/gnomAD/v4.1.1/genomes/genomes.bb",
    "gnomad_sv": f"{GBDB}/gnomAD/v4/structuralVariants/gnomad.v4.1.sv.non_neuro_controls.sites.bb",
    "trexplorer": f"{GBDB}/strVar/trexplorer.bb",
    "gtex_dapg": f"{GBDB}/gtex/eQtl/gtexDapg.bb",
    "mpravardb": f"{GBDB}/mpra/mpravardb/mpravardb.bb",
}
SOURCE_NOTES = {
    "maf": (
        "curated: HPRC release 1 Cactus alignment of 90 human assemblies on hg38 "
        "(UCSC hprc90way, Liao et al. 2023)"
    ),
    "hprc_arr": (
        "curated: HPRC arrangements relative to hg38 "
        "(UCSC hprcArrV1: duplications, inversions, deletions, insertions, double)"
    ),
    "hprc_v21_sv": (
        "curated: HPRC release 2 v2.1 structural variants from 233 assemblies "
        "(UCSC hprc2v21Sv, minigraph-cactus)"
    ),
    "gnomad_snv": "curated: gnomAD v4.1.1 genomes, allele frequencies (AN per site, up to 152,000 alleles)",
    "gnomad_sv": "curated: gnomAD v4.1 structural variants, allele frequencies (AN per site)",
    "trexplorer": "curated: TRExplorer v2 tandem-repeat catalogue, HPRC and TenK10K allele-size histograms",
    "gtex_dapg": "experimental: GTEx v8 fine-mapped cis-eQTLs (DAP-G, UCSC gtexEqtlDapg)",
    "mpravardb": "experimental: MPRAVarDB allele pairs tested in massively parallel reporter assays",
}
CACHE = Path("data/knowledge/human_panel")
USER_AGENT = "GenomeOS/0.9 (human panel; MAF ranges one chromosome at a time)"
CHUNK_BYTES = 32 << 20  # one range request carries at most this much MAF text
FETCH_THREADS = 4
INDEX_WINDOW = 10_000_000  # the offset table is asked for in windows of this size
GC_BIN = 100  # the hg38 row's GC is kept per 100 bp
INDEX_MARGIN = 30_000  # regional reads ask for offsets this far beyond each region
MAX_BLOCK_BYTES = 64 << 20  # no single MAF block is this large; a larger span means a hole in the index

# -- assembly status codes in a block ------------------------------------------------------------
PRESENT = ord("s")  # an `s` line: bases, gaps or N at each column
DELETED = ord("C")  # an `e` line with status C: contiguous on both sides, the block is not there
REPLACED = frozenset(b"InT")  # `e` status I, n or T: other sequence stands where the block is
MISSING = frozenset(b"M.")  # `e` status M (missing data) or no line at all
SAME = ord("=")  # in a site's allele string: the assembly carries hg38's base
ABSENT = ord(".")  # in a site's allele string: the assembly has no `s` line in this block
GAP = ord("-")
N_BASE = ord("N")
BASES = b"ACGT"

_NONZERO = re.compile(rb"[^\x00]+")


def assembly_of(src: str) -> str:
    """`GCA_018472595.1.JAHBCB010000033.1` -> `GCA_018472595.1`; `hs1.chr21` -> `hs1`."""
    if src.startswith("GCA_") or src.startswith("GCF_"):
        parts = src.split(".")
        return f"{parts[0]}.{parts[1]}"
    return src.split(".", 1)[0]


# ================================================================================================
# Range requests over kept-alive connections
# ================================================================================================
REQUEST_TIMEOUT = 60  # seconds a single range request may take before it is retried
REQUEST_RETRIES = 4
_POOL = threading.local()


class KeepAliveSource:
    """Random access to a URL by (offset, size) over one persistent HTTPS connection per host and thread.

    The standard reader (`bigwig._Source`) opens a new TLS connection for every range; profiling the
    chr21 run showed each connection costing 0.5 s on a quiet link and 5 s on a busy one, and a
    regional read issues hundreds. Same interface: `read`, `requests`, `bytes_fetched`, `close`.
    """

    remote = True

    def __init__(self, url: str, timeout: int = REQUEST_TIMEOUT, retries: int = REQUEST_RETRIES):
        u = urllib.parse.urlsplit(url)
        self.scheme, self.host, self.path = u.scheme, u.netloc, u.path + (f"?{u.query}" if u.query else "")
        self.timeout = timeout
        self.retries = retries
        self.requests = 0
        self.bytes_fetched = 0
        self.connections = 0

    def _conn(self, fresh: bool = False) -> http.client.HTTPConnection:
        pool = getattr(_POOL, "conns", None)
        if pool is None:
            pool = _POOL.conns = {}
        key = (self.scheme, self.host)
        conn = pool.get(key)
        if conn is None or fresh:
            if conn is not None:
                conn.close()
            cls = http.client.HTTPSConnection if self.scheme == "https" else http.client.HTTPConnection
            conn = pool[key] = cls(self.host, timeout=self.timeout)
            self.connections += 1
        return conn

    def get(self, headers: dict[str, str]) -> tuple[int, bytes]:
        delay = 1.0
        for attempt in range(self.retries):
            try:
                conn = self._conn(fresh=attempt > 0)
                conn.request("GET", self.path, headers={"User-Agent": USER_AGENT, **headers})
                r = conn.getresponse()
                data = r.read()
                if r.status >= 500:
                    raise OSError(f"server error {r.status}")
                return r.status, data
            except (OSError, http.client.HTTPException):
                if attempt == self.retries - 1:
                    raise
                time.sleep(delay)
                delay *= 2
        raise OSError("unreachable")  # pragma: no cover

    def read(self, offset: int, size: int) -> bytes:
        self.requests += 1
        self.bytes_fetched += size
        status, data = self.get({"Range": f"bytes={offset}-{offset + size - 1}"})
        if status == 200 and len(data) > size:
            raise OSError("server ignored the range request")
        if status not in (200, 206):
            raise OSError(f"range request failed ({status})")
        return data

    def close(self) -> None:
        pass


def open_source(source: str | Path):
    if isinstance(source, str) and source.startswith(("http://", "https://")):
        return KeepAliveSource(source)
    return _Source(source)


# ================================================================================================
# The MAF, streamed and distilled
# ================================================================================================
@dataclass
class MafRow:
    kind: str  # s, e or i
    assembly: str
    start: int = 0
    size: int = 0
    seq: bytes = b""
    status: str = ""
    left: str = ""
    left_count: int = 0


def parse_block(text: bytes) -> list[MafRow]:
    """The `s`, `e` and `i` lines of one MAF paragraph."""
    rows: list[MafRow] = []
    for line in text.split(b"\n"):
        if not line or line[1:2] != b" ":
            continue
        k = line[:1]
        if k == b"s":
            p = line.split()
            rows.append(MafRow("s", assembly_of(p[1].decode()), int(p[2]), int(p[3]), p[6].upper()))
        elif k == b"e":
            p = line.split()
            rows.append(MafRow("e", assembly_of(p[1].decode()), int(p[2]), int(p[3]), status=p[6].decode()))
        elif k == b"i":
            p = line.split()
            rows.append(MafRow("i", assembly_of(p[1].decode()), left=p[2].decode(), left_count=int(p[3])))
    return rows


def split_blocks(chunk: bytes) -> list[bytes]:
    """A run of whole MAF paragraphs cut into one bytes object each (header lines dropped)."""
    out = []
    for para in chunk.split(b"\n\n"):
        i = para.find(b"a ")
        if i < 0:
            continue
        if i and para[i - 1 : i] != b"\n":
            j = para.find(b"\na ")
            if j < 0:
                continue
            i = j + 1
        out.append(para[i:])
    return out


class Distiller:
    """Accumulates blocks of one chromosome into the store; the MAF text is discarded as it goes."""

    def __init__(self, chrom: str, reference: str = REFERENCE):
        self.chrom = chrom
        self.reference = reference
        self.assemblies: list[str] = []
        self._index: dict[str, int] = {}
        self.blocks: list[tuple[int, int, bytes, int]] = []  # start, end, status per assembly, extra copies
        self.sites: list[tuple[int, int, bytes]] = []  # pos, hg38 base, allele per assembly
        self.insertions: dict[int, dict[int, int]] = {}  # anchor -> assembly -> inserted bases (aligned)
        self.unaligned_insertions: dict[int, dict[int, int]] = {}  # anchor -> assembly -> non-aligning bases
        self.gc: dict[int, list[int]] = {}  # bin -> [G+C, A+C+G+T]
        self.blocks_seen = 0

    def _idx(self, assembly: str) -> int:
        i = self._index.get(assembly)
        if i is None:
            i = self._index[assembly] = len(self.assemblies)
            self.assemblies.append(assembly)
        return i

    def add(self, rows: list[MafRow]) -> None:
        ref = next((r for r in rows if r.kind == "s" and r.assembly == self.reference), None)
        if ref is None:
            return
        self.blocks_seen += 1
        seq = ref.seq
        width = len(seq)
        gaps = [m.start() for m in re.finditer(rb"-", seq)] if b"-" in seq else []
        refint = int.from_bytes(seq, "big")
        status: dict[int, int] = {}
        extra = 0
        diffs: list[tuple[int, int, int]] = []
        for r in rows:
            if r.assembly == self.reference:
                continue
            idx = self._idx(r.assembly)
            if r.kind == "s":
                if status.get(idx) == PRESENT:
                    extra += 1  # a second copy of the block in the same assembly
                    continue
                status[idx] = PRESENT
                if r.seq == seq or len(r.seq) != width:
                    continue
                d = (int.from_bytes(r.seq, "big") ^ refint).to_bytes(width, "big")
                s = r.seq
                for m in _NONZERO.finditer(d):
                    for col in range(m.start(), m.end()):
                        diffs.append((col, idx, s[col]))
            elif r.kind == "e":
                if status.get(idx) != PRESENT:
                    status[idx] = ord(r.status[:1] or "M")
            elif r.kind == "i" and r.left == "I" and r.left_count > 0:
                slot = self.unaligned_insertions.setdefault(ref.start - 1, {})
                slot[idx] = max(slot.get(idx, 0), r.left_count)
        n = len(self.assemblies)
        st = bytearray(b"." * n)
        for i, v in status.items():
            st[i] = v
        template = bytearray(b"." * n)
        for i, v in status.items():
            if v == PRESENT:
                template[i] = SAME
        block_sites: dict[int, bytearray] = {}
        for col, idx, ch in diffs:
            before = bisect_left(gaps, col) if gaps else 0
            if seq[col] == GAP:
                anchor = ref.start + col - before - 1
                slot = self.insertions.setdefault(anchor, {})
                slot[idx] = slot.get(idx, 0) + 1
                continue
            if seq[col] == N_BASE:
                continue
            pos = ref.start + col - before
            a = block_sites.get(pos)
            if a is None:
                a = block_sites[pos] = bytearray(template)
            a[idx] = ch
        bare = seq.replace(b"-", b"") if gaps else seq
        for pos in sorted(block_sites):
            self.sites.append((pos, bare[pos - ref.start], bytes(block_sites[pos])))
        self.blocks.append((ref.start, ref.start + ref.size, bytes(st), extra))
        # GC of the hg38 row, per bin
        p = ref.start
        while p < ref.start + ref.size:
            b = p // GC_BIN
            hi = min((b + 1) * GC_BIN, ref.start + ref.size)
            piece = bare[p - ref.start : hi - ref.start]
            g = piece.count(b"G") + piece.count(b"C")
            acgt = g + piece.count(b"A") + piece.count(b"T")
            slot = self.gc.setdefault(b, [0, 0])
            slot[0] += g
            slot[1] += acgt
            p = hi

    # -- persistence -----------------------------------------------------------------------------
    def save(self, directory: Path, meta: dict[str, Any] | None = None) -> Path:
        directory.mkdir(parents=True, exist_ok=True)
        n = len(self.assemblies)

        def pad(b: bytes) -> str:
            return (b + b"." * (n - len(b))).decode("ascii")

        with gzip.open(directory / "blocks.tsv.gz", "wt") as fh:
            for s, e, st, extra in sorted(self.blocks):
                fh.write(f"{s}\t{e}\t{pad(st)}\t{extra}\n")
        with gzip.open(directory / "sites.tsv.gz", "wt") as fh:
            for pos, base, alleles in sorted(self.sites):
                fh.write(f"{pos}\t{chr(base)}\t{pad(alleles)}\n")
        with gzip.open(directory / "insertions.tsv.gz", "wt") as fh:
            for kind, table in (("aligned", self.insertions), ("unaligned", self.unaligned_insertions)):
                for anchor in sorted(table):
                    cells = ",".join(f"{i}:{v}" for i, v in sorted(table[anchor].items()))
                    fh.write(f"{anchor}\t{kind}\t{cells}\n")
        with gzip.open(directory / "gc.tsv.gz", "wt") as fh:
            for b in sorted(self.gc):
                fh.write(f"{b}\t{self.gc[b][0]}\t{self.gc[b][1]}\n")
        info = {
            "chrom": self.chrom,
            "reference": self.reference,
            "assemblies": self.assemblies,
            "blocks": len(self.blocks),
            "sites": len(self.sites),
            "gc_bin": GC_BIN,
            **(meta or {}),
        }
        (directory / "meta.json").write_text(json.dumps(info, indent=1))
        return directory


def maf_index(chrom: str, start: int = 0, end: int | None = None, api=_api) -> list[tuple[int, int, int]]:
    """UCSC's table of MAF block offsets for a range: (hg38 start, hg38 end, byte offset), by offset."""
    items: dict[int, tuple[int, int, int]] = {}
    stop = end if end is not None else chromosome_size(chrom, api)
    lo = start
    while lo < stop:
        hi = min(lo + INDEX_WINDOW, stop)
        d = api(
            {
                "genome": "hg38",
                "track": MAF_TRACK,
                "chrom": chrom,
                "start": lo,
                "end": hi,
                "maxItemsOutput": 1_000_000,
            }
        )
        if "error" in d:
            raise OSError(d["error"])
        for it in d.get(MAF_TRACK, []):
            items[it["offset"]] = (it["chromStart"], it["chromEnd"], it["offset"])
        lo = hi
        time.sleep(0.3)
    return [items[o] for o in sorted(items)]


_SIZES: dict[str, int] = {}


def chromosome_size(chrom: str, api=None) -> int:
    if not _SIZES:
        _SIZES.update(_list_chromosomes())
    return _SIZES[chrom]


def _list_chromosomes() -> dict[str, int]:
    import urllib.request

    req = urllib.request.Request(
        "https://api.genome.ucsc.edu/list/chromosomes?genome=hg38", headers={"User-Agent": USER_AGENT}
    )
    with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310
        return {k: int(v) for k, v in json.load(r)["chromosomes"].items()}


def remote_size(url: str, retries: int = REQUEST_RETRIES) -> int:
    """The file's length from a HEAD request, retried over a fresh connection when the link stalls."""
    src = KeepAliveSource(url)
    delay = 2.0
    for attempt in range(retries):
        try:
            conn = src._conn(fresh=attempt > 0)
            conn.request("HEAD", src.path, headers={"User-Agent": USER_AGENT})
            r = conn.getresponse()
            r.read()
            return int(r.headers["Content-Length"])
        except (OSError, http.client.HTTPException, TypeError):
            if attempt == retries - 1:
                raise
            time.sleep(delay)
            delay *= 2
    raise OSError("unreachable")  # pragma: no cover


def byte_runs(
    index: list[tuple[int, int, int]],
    wanted: list[tuple[int, int]] | None,
    file_size: int,
    limit: int = CHUNK_BYTES,
) -> list[tuple[int, int]]:
    """Byte ranges (offset, size) covering the index entries that overlap `wanted` (all when None).

    A block runs from its offset to the next entry's offset; consecutive blocks are joined while
    the request stays under `limit`, so each range holds whole paragraphs.
    """
    offsets = [o for _, _, o in index] + [file_size]
    wanted = sorted(wanted) if wanted else None
    wstarts = [s for s, _ in wanted] if wanted else []
    wends = [e for _, e in wanted] if wanted else []
    runs: list[tuple[int, int]] = []
    for k, (s, e, o) in enumerate(index):
        if wanted is not None:
            i = bisect_right(wends, s)
            if not (i < len(wstarts) and wstarts[i] < e):
                continue
        size = offsets[k + 1] - o
        if size > MAX_BLOCK_BYTES:
            raise ValueError(
                f"block at offset {o} would read {size} bytes: the index around it is incomplete"
            )
        if runs and runs[-1][0] + runs[-1][1] == o and runs[-1][1] + size <= limit:
            runs[-1] = (runs[-1][0], runs[-1][1] + size)
        else:
            runs.append((o, size))
    return runs


def distil_chromosome(
    chrom: str,
    wanted: list[tuple[int, int]] | None = None,
    url: str | None = None,
    cache: Path = CACHE,
    threads: int = FETCH_THREADS,
    progress=None,
    name: str | None = None,
) -> Path:
    """Stream a chromosome's MAF (or the blocks over `wanted`) and save the distilled store."""
    url = url or MAF_URL.format(chrom=chrom)
    t0 = time.time()
    if wanted:
        # one offset query per region, widened so the block after the last wanted one is known
        merged: list[list[int]] = []
        for s, e in sorted(wanted):
            if merged and s - INDEX_MARGIN <= merged[-1][1] + INDEX_MARGIN:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        found: dict[int, tuple[int, int, int]] = {}
        for s, e in merged:
            for it in maf_index(chrom, max(0, s - INDEX_MARGIN), e + INDEX_MARGIN):
                found[it[2]] = it
        index = [found[o] for o in sorted(found)]
    else:
        index = maf_index(chrom)
    size = remote_size(url)
    runs = byte_runs(index, wanted, size)
    total = sum(n for _, n in runs)
    dist = Distiller(chrom)
    done = 0

    def fetch(run: tuple[int, int]) -> bytes:
        return KeepAliveSource(url, timeout=300).read(run[0], run[1])

    with ThreadPoolExecutor(max_workers=threads) as pool:
        for i in range(0, len(runs), threads * 2):
            batch = runs[i : i + threads * 2]
            for run, chunk in zip(batch, pool.map(fetch, batch), strict=True):
                for para in split_blocks(chunk):
                    dist.add(parse_block(para))
                done += run[1]
                if progress:
                    progress(done, total, dist.blocks_seen, time.time() - t0)
    meta = {
        "source": url,
        "index_entries": len(index),
        "mb_fetched": round(total / 1e6, 1),
        "requests": len(runs),
        "seconds": round(time.time() - t0, 1),
        "wanted_intervals": len(wanted) if wanted else None,
    }
    return dist.save(cache / (name or chrom), meta)


# ================================================================================================
# The store, read back
# ================================================================================================
@dataclass
class SiteCounts:
    pos: int
    ref: int
    present: int  # assemblies with an `s` line here
    called: int  # of those, a base other than N
    gap: int
    alt: dict[int, int]  # alternative base -> assemblies
    minor: int  # assemblies not carrying the commonest state (bases and gap)


class Panel:
    """One chromosome's distilled store: blocks, variable columns, insertions, GC."""

    def __init__(self, directory: Path):
        self.directory = Path(directory)
        self.meta = json.loads((self.directory / "meta.json").read_text())
        self.chrom = self.meta["chrom"]
        self.assemblies: list[str] = self.meta["assemblies"]
        self.n = len(self.assemblies)
        self.block_start: list[int] = []
        self.block_end: list[int] = []
        self.block_status: list[bytes] = []
        self.block_extra: list[int] = []
        with gzip.open(self.directory / "blocks.tsv.gz", "rt") as fh:
            for line in fh:
                s, e, st, extra = line.rstrip("\n").split("\t")
                self.block_start.append(int(s))
                self.block_end.append(int(e))
                self.block_status.append(st.encode("ascii"))
                self.block_extra.append(int(extra))
        self.site_pos: list[int] = []
        self.site_ref: list[int] = []
        self.site_alleles: list[bytes] = []
        with gzip.open(self.directory / "sites.tsv.gz", "rt") as fh:
            for line in fh:
                p, b, a = line.rstrip("\n").split("\t")
                self.site_pos.append(int(p))
                self.site_ref.append(ord(b))
                self.site_alleles.append(a.encode("ascii"))
        self.ins_aligned: dict[int, dict[int, int]] = {}
        self.ins_unaligned: dict[int, dict[int, int]] = {}
        with gzip.open(self.directory / "insertions.tsv.gz", "rt") as fh:
            for line in fh:
                anchor, kind, cells = line.rstrip("\n").split("\t")
                table = self.ins_aligned if kind == "aligned" else self.ins_unaligned
                table[int(anchor)] = {int(i): int(v) for i, v in (c.split(":") for c in cells.split(","))}
        self.ins_anchors = sorted(set(self.ins_aligned) | set(self.ins_unaligned))
        self.gc: dict[int, tuple[int, int]] = {}
        with gzip.open(self.directory / "gc.tsv.gz", "rt") as fh:
            for line in fh:
                b, g, t = line.split("\t")
                self.gc[int(b)] = (int(g), int(t))
        self._counts: dict[int, SiteCounts] = {}

    # -- elementary reads ------------------------------------------------------------------------
    def blocks_over(self, start: int, end: int) -> range:
        i = bisect_right(self.block_end, start)
        j = bisect_left(self.block_start, end)
        return range(i, max(i, j))

    def sites_over(self, start: int, end: int) -> range:
        return range(bisect_left(self.site_pos, start), bisect_left(self.site_pos, end))

    def counts(self, k: int) -> SiteCounts:
        c = self._counts.get(k)
        if c is not None:
            return c
        a = self.site_alleles[k]
        ref = self.site_ref[k]
        absent = a.count(ABSENT)
        present = self.n - absent
        gap = a.count(GAP)
        nb = a.count(N_BASE)
        alt = {b: a.count(b) for b in BASES if b != ref and a.count(b)}
        states = [present - gap - nb - sum(alt.values()), gap, *alt.values()]
        called = present - nb
        minor = called - max(states)
        c = SiteCounts(self.site_pos[k], ref, present, called - gap, gap, alt, minor)
        self._counts[k] = c
        return c

    def gc_fraction(self, intervals: list[tuple[int, int]]) -> float | None:
        g = t = 0
        for s, e in intervals:
            for b in range(s // GC_BIN, (e - 1) // GC_BIN + 1):
                v = self.gc.get(b)
                if v:
                    g += v[0]
                    t += v[1]
        return g / t if t else None

    def aligned_bases(self, start: int, end: int) -> int:
        return sum(
            min(end, self.block_end[i]) - max(start, self.block_start[i])
            for i in self.blocks_over(start, end)
        )


# ================================================================================================
# bigBed over HTTP ranges
# ================================================================================================
BIGBED_MAGIC = 0x8789F2EB


class BigBed(BigWig):
    """A bigBed opened by URL or path; items of one chromosome range come back as dicts.

    The chromosome tree, the R-tree walk and the coalesced range requests are the bigWig reader's
    (attribution/bigwig.py); only the header magic, the autoSql field names and the record format
    differ.
    """

    def __init__(self, source: str | Path):  # noqa: D107 - the parent's reader with bigBed's header
        self.src = open_source(source)
        h = self.src.read(0, 64)
        magic, version, zoom, ct, fd, fi, fc, dfc, asql, _ts, ubs, _res = struct.unpack("<IHHQQQHHQQIQ", h)
        if magic != BIGBED_MAGIC:
            raise ValueError(f"not a little-endian bigBed: magic {magic:#x}")
        self.header = Header(version, zoom, ct, fd, fi, ubs)
        self.field_count = fc
        self.fields = self._autosql(asql) if asql else []
        self.chroms: dict[str, tuple[int, int]] = {}
        self._read_chrom_tree()
        (_m, self._rtree_block_size, _n, _sc, _sb, _ec, _eb, _eo, self._items_per_slot, _r) = struct.unpack(
            "<IIQIIIIQII", self.src.read(fi, 48)
        )
        self._node_cache: dict[int, bytes] = {}

    def _autosql(self, offset: int) -> list[str]:
        raw = self.src.read(offset, 16_000)
        text = raw.split(b"\0", 1)[0].decode("ascii", "replace")
        names = []
        for line in text.splitlines():
            m = re.match(r"\s*[\w\[\]\s]+?\s+(\w+)\s*;", line)
            if m and not line.strip().startswith(("table", "(", ")", '"')):
                names.append(m.group(1))
        return names

    def query(self, chrom: str, intervals: list[tuple[int, int]], limit_bytes: int = 16 << 20) -> list[dict]:
        """Every item overlapping any of the sorted, non-overlapping intervals, once."""
        if not intervals or chrom not in self.chroms:
            return []
        ivs = sorted(intervals)
        cid = self.chroms[chrom][0]
        starts = [s for s, _ in ivs]
        ends = [e for _, e in ivs]
        out: list[dict] = []
        seen: set[tuple[int, int, str]] = set()
        for g in coalesce(self.leaf_items(chrom, ivs), max_size=limit_bytes):
            base = g[0].offset
            blob = self.src.read(base, g[-1].offset + g[-1].size - base)
            for it in g:
                raw = blob[it.offset - base : it.offset - base + it.size]
                if self.header.uncompress_buf_size:
                    raw = zlib.decompress(raw)
                p = 0
                while p + 12 <= len(raw):
                    c, s, e = struct.unpack_from("<III", raw, p)
                    z = raw.index(b"\0", p + 12)
                    rest = raw[p + 12 : z].decode("utf-8", "replace")
                    p = z + 1
                    if c != cid:
                        continue
                    i = bisect_right(ends, s)
                    if not (i < len(starts) and starts[i] < e):
                        continue
                    key = (s, e, rest[:64])
                    if key in seen:
                        continue
                    seen.add(key)
                    vals = rest.split("\t") if rest else []
                    row: dict[str, Any] = {"chrom": chrom, "chromStart": s, "chromEnd": e}
                    for name, v in zip(self.fields[3:], vals, strict=False):
                        row[name] = v
                    out.append(row)
        return out


_OPEN_BIGBEDS: dict[str, BigBed] = {}


def open_bigbed(key: str) -> BigBed:
    """One reader per track for the life of the process: header, autoSql and chromosome tree once."""
    bb = _OPEN_BIGBEDS.get(key)
    if bb is None:
        bb = _OPEN_BIGBEDS[key] = BigBed(BIGBEDS[key])
    return bb


# ================================================================================================
# Events: what varies inside an interval, one entry per variant, with every assembly's state
# ================================================================================================
MIN_RECURRING = 2  # a state carried by two or more assemblies recurs; one is a singleton or an error


@dataclass
class Event:
    kind: str  # snv, deletion, insertion, structural
    start: int
    end: int
    tokens: list[Any]  # per assembly: a state, or None where the assembly does not inform
    ref_token: Any

    @property
    def informative(self) -> int:
        return sum(t is not None for t in self.tokens)

    @property
    def minor(self) -> int:
        """Assemblies not carrying the commonest state; an absent assembly ('*') belongs to the
        structural event that removed it and is not counted again here."""
        c = Counter(t for t in self.tokens if t is not None and t != "*")
        return sum(c.values()) - max(c.values()) if c else 0

    def states(self) -> dict[str, int]:
        c = Counter(str(t) for t in self.tokens if t is not None)
        return dict(c.most_common())


def _block_at(panel: Panel, pos: int) -> int | None:
    i = bisect_right(panel.block_start, pos) - 1
    if i >= 0 and panel.block_end[i] > pos:
        return i
    return None


def _absent_token(code: int) -> str | None:
    if code == DELETED:
        return "*"
    if code in REPLACED:
        return "*"
    return None


def events(panel: Panel, intervals: list[tuple[int, int]]) -> list[Event]:
    """Every variable column, deletion run, insertion and absent block inside the intervals."""
    n = panel.n
    out: list[Event] = []
    for s, e in sorted(intervals):
        # blocks where some assembly is deleted or replaced: one structural event per run
        prev: Event | None = None
        for bi in panel.blocks_over(s, e):
            st = panel.block_status[bi]
            if st.count(PRESENT) == n:
                prev = None
                continue
            toks: list[Any] = []
            for code in st:
                if code == PRESENT:
                    toks.append("=")
                elif code == DELETED:
                    toks.append("D")
                elif code in REPLACED:
                    toks.append("R")
                else:
                    toks.append(None)
            lo, hi = max(s, panel.block_start[bi]), min(e, panel.block_end[bi])
            if prev is not None and prev.end == lo and prev.tokens == toks:
                prev.end = hi
                continue
            if any(t in ("D", "R") for t in toks):
                prev = Event("structural", lo, hi, toks, "=")
                out.append(prev)
            else:
                prev = None
        # variable columns; runs of gaps shared by the same assemblies are one deletion
        run: Event | None = None
        for k in panel.sites_over(s, e):
            pos = panel.site_pos[k]
            a = panel.site_alleles[k]
            bi = _block_at(panel, pos)
            st = panel.block_status[bi] if bi is not None else b"." * n
            toks = []
            has_alt = has_gap = False
            for i, c in enumerate(a):
                if c == SAME:
                    toks.append("=")
                elif c == GAP:
                    toks.append("-")
                    has_gap = True
                elif c == ABSENT:
                    toks.append(_absent_token(st[i]))
                elif c == N_BASE:
                    toks.append(None)
                else:
                    toks.append(chr(c))
                    has_alt = True
            if has_gap and not has_alt:
                if run is not None and run.end == pos and run.tokens == toks:
                    run.end = pos + 1
                    continue
                run = Event("deletion", pos, pos + 1, toks, "=")
                out.append(run)
                continue
            run = None
            if not has_alt and not has_gap:
                continue
            out.append(Event("snv", pos, pos + 1, toks, "="))
        # insertions after an hg38 base inside the intervals
        for j in range(bisect_left(panel.ins_anchors, s), bisect_left(panel.ins_anchors, e)):
            anchor = panel.ins_anchors[j]
            lengths: dict[int, int] = {}
            for table in (panel.ins_aligned, panel.ins_unaligned):
                for i, v in table.get(anchor, {}).items():
                    lengths[i] = lengths.get(i, 0) + v
            bi = _block_at(panel, anchor)
            st = panel.block_status[bi] if bi is not None else b"." * n
            toks = []
            for i in range(n):
                if st[i] == PRESENT:
                    toks.append(lengths.get(i, 0))
                elif i in lengths:
                    toks.append(lengths[i])
                else:
                    toks.append(_absent_token(st[i]))
            out.append(Event("insertion", anchor, anchor + 1, toks, 0))
    out.sort(key=lambda ev: (ev.start, ev.kind))
    return out


# ================================================================================================
# Blocks: presence, identity, structure, diversity
# ================================================================================================
def heterozygosity(tokens: list[Any]) -> float:
    """Expected heterozygosity among the base states of one column (gaps and absences left out)."""
    c = Counter(t for t in tokens if t is not None and t not in ("-", "*"))
    m = sum(c.values())
    if m < 2:
        return 0.0
    return m / (m - 1) * (1.0 - sum((v / m) ** 2 for v in c.values()))


def measure(panel: Panel, intervals: list[tuple[int, int]], evs: list[Event] | None = None) -> dict[str, Any]:
    """The panel's numbers over a set of hg38 intervals (one block, or one gene's coding exons)."""
    n = panel.n
    ivs = sorted(intervals)
    bases = sum(e - s for s, e in ivs)
    aligned = 0
    cells = Counter()
    per_assembly = [0] * n
    for s, e in ivs:
        for bi in panel.blocks_over(s, e):
            ov = min(e, panel.block_end[bi]) - max(s, panel.block_start[bi])
            if ov <= 0:
                continue
            aligned += ov
            st = panel.block_status[bi]
            cells["present"] += ov * st.count(PRESENT)
            cells["deleted"] += ov * st.count(DELETED)
            cells["replaced"] += ov * sum(st.count(c) for c in REPLACED)
            for i, code in enumerate(st):
                if code == PRESENT:
                    per_assembly[i] += ov
    gaps = ns = alt = 0
    for s, e in ivs:
        for k in panel.sites_over(s, e):
            a = panel.site_alleles[k]
            g = a.count(GAP)
            nb = a.count(N_BASE)
            gaps += g
            ns += nb
            alt += sum(a.count(b) for b in BASES)
            if g:
                bi = _block_at(panel, panel.site_pos[k])
                if bi is not None:
                    for i, c in enumerate(a):
                        if c == GAP:
                            per_assembly[i] -= 1
    evs = events(panel, ivs) if evs is None else evs
    informative_cells = cells["present"] + cells["deleted"] + cells["replaced"]
    called = cells["present"] - gaps - ns
    by_kind: dict[str, dict[str, int]] = {}
    pi_sum = 0.0
    for ev in evs:
        k = by_kind.setdefault(ev.kind, {"events": 0, "recurring": 0})
        k["events"] += 1
        if ev.minor >= MIN_RECURRING:
            k["recurring"] += 1
        if ev.kind == "snv":
            pi_sum += heterozygosity(ev.tokens)
    touched, touched_recurring = clip_touched(evs, ivs)
    recurring = sum(v["recurring"] for v in by_kind.values())
    return {
        "bases": bases,
        "aligned_bases": aligned,
        "aligned_share": round(aligned / bases, 4) if bases else None,
        "presence": round((cells["present"] - gaps) / informative_cells, 4) if informative_cells else None,
        "deleted_share": round((cells["deleted"] + gaps) / informative_cells, 5)
        if informative_cells
        else None,
        "replaced_share": round(cells["replaced"] / informative_cells, 5) if informative_cells else None,
        "missing_share": round(1 - informative_cells / (aligned * n), 4) if aligned else None,
        "assemblies_aligned_any": sum(v > 0 for v in per_assembly),
        "assemblies_aligned_90": sum(v >= 0.9 * aligned for v in per_assembly) if aligned else 0,
        "identity": round(1 - alt / called, 6) if called else None,
        "events": by_kind,
        "recurring_events": recurring,
        "recurring_per_kb": round(1000 * recurring / aligned, 3) if aligned else None,
        "pi": round(pi_sum / aligned, 6) if aligned else None,
        "touched_share": round(touched / bases, 4) if bases else None,
        "touched_recurring_share": round(touched_recurring / bases, 4) if bases else None,
    }


# ================================================================================================
# Units: the value domain
# ================================================================================================
UNIT_WIDTH = 200  # bases per unit when blocks are tiled
UNIT_MIN_ALIGNED = 0.5  # the panel must align at least half a unit's bases
UNIT_MIN_INFORMATIVE = 0.8  # and at least 80% of the assemblies must inform every recurring event
FIXED_MIN_SHARE = 0.95  # one value in 95% of the assemblies: fixed
STORAGE_MIN_COVER = 0.8  # recurring values carried by 80% or more of the assemblies
STORAGE_MAX_VALUES = 8  # and no more than this many recurring values: a column with a small domain


def entropy_bits(counts: list[int]) -> float:
    m = sum(counts)
    return -sum(c / m * math.log2(c / m) for c in counts if c) if m else 0.0


def domain(panel: Panel, intervals: list[tuple[int, int]], evs: list[Event] | None = None) -> dict[str, Any]:
    """Each assembly's value over the unit's recurring events, and the distribution of values."""
    ivs = sorted(intervals)
    bases = sum(e - s for s, e in ivs)
    aligned = sum(panel.aligned_bases(s, e) for s, e in ivs)
    evs = events(panel, ivs) if evs is None else evs
    rec = [ev for ev in evs if ev.minor >= MIN_RECURRING]
    values: list[tuple | None] = []
    for i in range(panel.n):
        toks = tuple(ev.tokens[i] for ev in rec)
        values.append(None if any(t is None for t in toks) else toks)
    inf = [v for v in values if v is not None]
    counts = Counter(inf)
    ordered = counts.most_common()
    m = len(inf)
    recurring = [(v, c) for v, c in ordered if c >= MIN_RECURRING]
    out: dict[str, Any] = {
        "bases": bases,
        "aligned_share": round(aligned / bases, 3) if bases else 0.0,
        "informative": m,
        "recurring_events": len(rec),
        "values": len(counts),
        "recurring_values": len(recurring),
        "top_share": round(ordered[0][1] / m, 3) if m else None,
        "recurring_cover": round(sum(c for _, c in recurring) / m, 3) if m else None,
        "entropy_bits": round(entropy_bits([c for _, c in ordered]), 3),
        "effective_values": round(1 / sum((c / m) ** 2 for _, c in ordered), 2) if m else None,
    }
    out["class"] = unit_class(out, panel.n)
    out["_events"] = rec
    out["_values"] = values
    return out


def unit_class(
    d: dict[str, Any],
    n: int,
    fixed_min: float = FIXED_MIN_SHARE,
    cover_min: float = STORAGE_MIN_COVER,
    max_values: int = STORAGE_MAX_VALUES,
) -> str:
    if d["aligned_share"] < UNIT_MIN_ALIGNED or d["informative"] < UNIT_MIN_INFORMATIVE * n:
        return "unplaced"
    if d["top_share"] >= fixed_min:
        return "fixed"
    if d["recurring_cover"] >= cover_min and d["recurring_values"] <= max_values:
        return "storage"
    return "hypervariable"


def describe_values(d: dict[str, Any], limit: int = 8) -> list[dict[str, Any]]:
    """The commonest values of a unit, each as its departures from hg38 and how many carry it."""
    rec: list[Event] = d["_events"]
    counts = Counter(v for v in d["_values"] if v is not None)
    m = sum(counts.values())
    out = []
    for v, c in counts.most_common(limit):
        diffs = []
        for ev, t in zip(rec, v, strict=True):
            if t == ev.ref_token:
                continue
            if ev.kind == "snv":
                diffs.append(
                    f"{ev.start + 1}{'del' if t == '-' else '>' + str(t)}"
                    if t != "*"
                    else f"{ev.start + 1}:absent"
                )
            elif ev.kind == "deletion":
                diffs.append(
                    f"{ev.start + 1}-{ev.end}del{ev.end - ev.start}" if t == "-" else f"{ev.start + 1}:absent"
                )
            elif ev.kind == "insertion":
                diffs.append(f"{ev.start + 1}ins{t}" if t != "*" else f"{ev.start + 1}:absent")
            else:
                diffs.append(f"{ev.start + 1}-{ev.end}{'deleted' if t == 'D' else 'replaced'}")
        out.append({"assemblies": c, "share": round(c / m, 3), "differs_from_hg38": diffs or ["hg38"]})
    return out


def alternatives(d: dict[str, Any]) -> list[dict[str, Any]]:
    """The recurring alternatives of a unit, one per event: substitution, indel length or absence."""
    out = []
    for ev in d["_events"]:
        row: dict[str, Any] = {"kind": ev.kind, "start": ev.start, "end": ev.end, "states": ev.states()}
        if ev.kind == "insertion":
            lengths = [t for t in ev.tokens if isinstance(t, int)]
            row["inserted_lengths"] = dict(Counter(lengths).most_common(12))
        out.append(row)
    return out


def length_domain(panel: Panel, intervals: list[tuple[int, int]], evs: list[Event]) -> dict[int, int]:
    """Net bases gained or lost relative to hg38 per assembly over the unit: the length allele."""
    net = [0] * panel.n
    ok = [True] * panel.n
    for ev in evs:
        for i, t in enumerate(ev.tokens):
            if t is None:
                ok[i] = False
            elif ev.kind == "insertion" and isinstance(t, int):
                net[i] += t
            elif (ev.kind == "deletion" and t == "-") or (ev.kind == "structural" and t in ("D", "R")):
                net[i] -= ev.end - ev.start
    return dict(Counter(v for v, good in zip(net, ok, strict=True) if good).most_common())


def length_class(
    lengths: dict[int, int],
    n: int,
    fixed_min: float = FIXED_MIN_SHARE,
    cover_min: float = STORAGE_MIN_COVER,
    max_values: int = STORAGE_MAX_VALUES,
) -> dict[str, Any]:
    """The same rule on the length alleles alone: a repeat can be a column of copy numbers while the
    sequence of every copy differs."""
    m = sum(lengths.values())
    counts = sorted(lengths.values(), reverse=True)
    recurring = [c for c in counts if c >= MIN_RECURRING]
    out: dict[str, Any] = {
        "informative": m,
        "values": len(counts),
        "recurring_values": len(recurring),
        "top_share": round(counts[0] / m, 3) if m else None,
        "recurring_cover": round(sum(recurring) / m, 3) if m else None,
        "effective_values": round(1 / sum((c / m) ** 2 for c in counts), 2) if m else None,
    }
    if m < UNIT_MIN_INFORMATIVE * n:
        out["class"] = "unplaced"
    elif out["top_share"] >= fixed_min:
        out["class"] = "fixed"
    elif out["recurring_cover"] >= cover_min and len(recurring) <= max_values:
        out["class"] = "storage"
    else:
        out["class"] = "hypervariable"
    return out


def tile(intervals: list[tuple[int, int]], width: int = UNIT_WIDTH) -> list[list[tuple[int, int]]]:
    """Cut a unit's bases, in order, into pieces of `width` bases; a short last piece joins the one before."""
    pieces: list[list[tuple[int, int]]] = []
    cur: list[tuple[int, int]] = []
    filled = 0
    for s, e in sorted(intervals):
        p = s
        while p < e:
            take = min(e - p, width - filled)
            cur.append((p, p + take))
            filled += take
            p += take
            if filled == width:
                pieces.append(cur)
                cur, filled = [], 0
    if cur:
        if pieces and filled < width // 2:
            pieces[-1].extend(cur)
        else:
            pieces.append(cur)
    return pieces


# ================================================================================================
# Background: recurring events per aligned base, by GC, from the chromosome itself
# ================================================================================================
BACKGROUND_BIN = 1000  # the background is counted per kilobase of the hg38 grid
GC_STRATA = (0.35, 0.40, 0.45, 0.50, 0.55, 0.60)  # stratum edges for GC matching
CORE_MAX_RATIO = 0.5  # recurring events at half the GC-matched background or less: depleted
CORE_MAX_P = 0.01  # and a Poisson lower tail at or below this
MIN_EXPECTED = 5.0  # fewer expected events than this and a block is too short to call
PRESENT_MIN = 0.95  # presence below this (deleted or replaced in 5% of assembly-bases): polymorphic
LINEAGE_MAX = 0.5  # presence below this: most assemblies lack the sequence
STRUCTURAL_MAX = 0.10  # a tenth of the bases touched by a recurring deletion, replacement or insertion
ALIGNED_MIN = 0.5  # the panel aligns less than half the block: unplaced
MISSING_MAX = 0.25  # or a quarter of the assembly-bases are missing data


def gc_stratum(gc: float | None) -> int | None:
    if gc is None:
        return None
    return bisect_right(GC_STRATA, gc)


def poisson_tails(observed: int, expected: float) -> tuple[float, float]:
    """P(X <= observed) and P(X >= observed) for a Poisson with the expected mean."""
    if expected <= 0:
        return 1.0, 1.0 if observed == 0 else 0.0
    if expected > 500:
        z = (observed + 0.5 - expected) / math.sqrt(expected)
        lo = 0.5 * math.erfc(-z / math.sqrt(2))
        z2 = (observed - 0.5 - expected) / math.sqrt(expected)
        hi = 0.5 * math.erfc(z2 / math.sqrt(2))
        return lo, hi
    logp = -expected
    cdf = 0.0
    term = math.exp(logp)
    for k in range(observed + 1):
        if k:
            term *= expected / k
        cdf += term
    below = min(1.0, cdf)
    above = min(1.0, 1.0 - cdf + term)
    return below, above


@dataclass
class Background:
    """Recurring events and aligned bases per kilobase bin, and the density per stratum.

    A stratum is the GC stratum of the kilobase, and with replication timing given, the pair of
    GC and timing strata, so an expectation can be matched on either.
    """

    bins: dict[int, list[Any]] = field(default_factory=dict)  # bin -> [aligned, recurring, gc]
    density: dict[Any, float] = field(default_factory=dict)
    overall: float = 0.0
    rt: dict[int, float] | None = None
    edges: list[float] | None = None

    def stratum(self, b: int) -> Any:
        row = self.bins.get(b)
        gs = gc_stratum(row[2]) if row else None
        if self.rt is None:
            return gs
        return (gs, rt_stratum(self.rt.get(b), self.edges or []))

    def rate(self, b: int) -> float:
        return self.density.get(self.stratum(b), self.overall)

    def expected(self, panel: Panel, intervals: list[tuple[int, int]]) -> float:
        total = 0.0
        for s, e in intervals:
            for b in range(s // BACKGROUND_BIN, (e - 1) // BACKGROUND_BIN + 1):
                lo, hi = max(s, b * BACKGROUND_BIN), min(e, (b + 1) * BACKGROUND_BIN)
                if lo < hi:
                    total += panel.aligned_bases(lo, hi) * self.rate(b)
        return total


def build_background(
    panel: Panel,
    exclude: list[tuple[int, int]],
    evs: list[Event],
    rt: dict[int, float] | None = None,
    edges: list[float] | None = None,
) -> Background:
    """The chromosome's own rate of recurring events per aligned base by stratum, outside `exclude`."""
    ex = merge_intervals(exclude)
    ex_starts = [s for s, _ in ex]
    ex_ends = [e for _, e in ex]

    def excluded(pos: int) -> bool:
        i = bisect_right(ex_ends, pos)
        return i < len(ex_starts) and ex_starts[i] <= pos

    bg = Background(rt=rt, edges=edges)
    for bi in range(len(panel.block_start)):
        s, e = panel.block_start[bi], panel.block_end[bi]
        for b in range(s // BACKGROUND_BIN, (e - 1) // BACKGROUND_BIN + 1):
            lo, hi = max(s, b * BACKGROUND_BIN), min(e, (b + 1) * BACKGROUND_BIN)
            cut = 0
            j = bisect_right(ex_ends, lo)
            while j < len(ex) and ex_starts[j] < hi:
                cut += max(0, min(hi, ex_ends[j]) - max(lo, ex_starts[j]))
                j += 1
            row = bg.bins.setdefault(b, [0, 0, None])
            row[0] += hi - lo - cut
    for b, row in bg.bins.items():
        row[2] = panel.gc_fraction([(b * BACKGROUND_BIN, (b + 1) * BACKGROUND_BIN)])
    for ev in evs:
        if ev.minor < MIN_RECURRING or excluded(ev.start):
            continue
        row = bg.bins.get(ev.start // BACKGROUND_BIN)
        if row is not None:
            row[1] += 1
    num: Counter = Counter()
    den: Counter = Counter()
    for b, row in bg.bins.items():
        if row[0] < BACKGROUND_BIN / 2 or row[2] is None:
            continue
        st = bg.stratum(b)
        num[st] += row[1]
        den[st] += row[0]
    bg.density = {st: num[st] / den[st] for st in den if den[st]}
    bg.overall = sum(num.values()) / sum(den.values()) if den else 0.0
    return bg


def block_class(m: dict[str, Any], observed: int, expected: float) -> dict[str, Any]:
    """Core, variable, polymorphic, lineage-restricted or unplaced, with the evidence and a confidence."""
    ratio = observed / expected if expected else None
    below, above = poisson_tails(observed, expected)
    ev = {
        "recurring_events": observed,
        "expected": round(expected, 1),
        "ratio": round(ratio, 3) if ratio is not None else None,
        "p_depleted": round(below, 5),
        "p_enriched": round(above, 5),
    }
    presence = m.get("presence")
    if (
        (m.get("aligned_share") or 0) < ALIGNED_MIN
        or (m.get("missing_share") or 0) > MISSING_MAX
        or presence is None
    ):
        return {"class": "unplaced", "reason": "not aligned", "confidence": 0.0, **ev}
    if presence < LINEAGE_MAX:
        conf = 0.3 + 0.3 * min(1.0, (LINEAGE_MAX - presence) / LINEAGE_MAX)
        return {"class": "lineage_restricted", "confidence": round(conf, 2), **ev}
    if presence < PRESENT_MIN or (m.get("touched_recurring_share") or 0) >= STRUCTURAL_MAX:
        gap = max(PRESENT_MIN - presence, (m.get("touched_recurring_share") or 0) - STRUCTURAL_MAX)
        return {"class": "polymorphic", "confidence": round(0.3 + min(0.3, gap * 3), 2), **ev}
    if expected < MIN_EXPECTED:
        return {"class": "unplaced", "reason": "too short to call", "confidence": 0.0, **ev}
    if ratio is not None and ratio <= CORE_MAX_RATIO and below <= CORE_MAX_P:
        conf = 0.3 + 0.3 * min(1.0, (CORE_MAX_RATIO - ratio) / CORE_MAX_RATIO + (1 if below < 1e-4 else 0))
        return {"class": "core", "confidence": round(min(conf, 0.6), 2), **ev}
    conf = 0.3 + 0.3 * min(1.0, abs((ratio or 0) - CORE_MAX_RATIO) / CORE_MAX_RATIO)
    return {"class": "variable", "confidence": round(min(conf, 0.6), 2), **ev}


# ================================================================================================
# Replication timing: ENCODE Repli-seq, lifted from hg19 per kilobase
# ================================================================================================
REPLI_BASE = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq"
REPLI_CELLS = {
    "BG02ES": "Bg02es",
    "BJ": "Bj",
    "GM12878": "Gm12878",
    "HeLa-S3": "Helas3",
    "HepG2": "Hepg2",
    "HUVEC": "Huvec",
    "IMR-90": "Imr90",
    "K562": "K562",
    "MCF-7": "Mcf7",
    "NHEK": "Nhek",
    "SK-N-SH": "Sknsh",
}
LIFT_CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
RT_STRATA = 3  # early, middle and late thirds of the chromosome's kilobases
RT_EVIDENCE = (
    "measured: ENCODE UW Repli-seq wavelet-smoothed signal (hg19, 11 cell lines, mean; high = early), "
    "lifted to hg38 per kilobase through UCSC's hg38ToHg19 chain"
)


def repli_url(cell: str) -> str:
    return f"{REPLI_BASE}/wgEncodeUwRepliSeq{REPLI_CELLS[cell]}WaveSignalRep1.bigWig"


@functools.lru_cache(maxsize=32)
def load_chain(chrom: str, path: Path | None = None) -> list[tuple[int, int, str, int, bool, int]]:
    """Ungapped blocks of the hg38->hg19 chain for one hg38 chromosome: (t_start, t_end, q_chrom,
    q_start, q_minus, q_size), sorted by t_start. The chain file is fetched once into the cache."""
    import urllib.request

    path = path or CACHE / "hg38ToHg19.over.chain.gz"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(LIFT_CHAIN_URL, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=300) as r:  # noqa: S310
            path.write_bytes(r.read())
    out = []
    keep = False
    t = q = 0
    q_chrom, q_minus, q_size = "", False, 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("chain"):
                p = line.split()
                keep = p[2] == chrom
                t, q = int(p[5]), int(p[10])
                q_chrom, q_size, q_minus = p[7], int(p[8]), p[9] == "-"
                continue
            if not keep or not line.strip():
                continue
            p = line.split()
            size = int(p[0])
            out.append((t, t + size, q_chrom, q, q_minus, q_size))
            if len(p) == 3:
                t += size + int(p[1])
                q += size + int(p[2])
    out.sort()
    return out


def lift_point(
    chain: list[tuple[int, int, str, int, bool, int]], starts: list[int], pos: int
) -> tuple[str, int] | None:
    i = bisect_right(starts, pos) - 1
    if i < 0:
        return None
    t0, t1, qc, q0, minus, qsize = chain[i]
    if pos >= t1:
        return None
    q = q0 + (pos - t0)
    return (qc, qsize - q - 1) if minus else (qc, q)


TIMING_CACHE = CACHE / "timing"


def repli_local(cell: str, cache: Path = CACHE) -> Path:
    """The cell's Repli-seq bigWig (about 10 MB) fetched once whole: a kilobase-resolution track is
    cheaper to read locally than by hundreds of ranges."""
    d = cache / "repliseq"
    d.mkdir(parents=True, exist_ok=True)
    p = d / Path(urllib.parse.urlsplit(repli_url(cell)).path).name
    if not p.exists():
        status, data = KeepAliveSource(repli_url(cell), timeout=300).get({})
        if status != 200:
            raise OSError(f"Repli-seq download failed ({status})")
        tmp = p.with_suffix(".part")
        tmp.write_bytes(data)
        tmp.rename(p)
    return p


def replication_timing(
    chrom: str, bins: list[int], cells: list[str] | None = None, chain=None, progress=None
) -> tuple[dict[int, float], dict[str, Any]]:
    """Mean Repli-seq signal per hg38 kilobase bin (centre lifted to hg19), and what it cost."""
    by, cost = replication_timing_many(
        {chrom: bins}, cells, {chrom: chain} if chain is not None else None, progress
    )
    return by.get(chrom, {}), cost


def replication_timing_many(
    wanted: dict[str, list[int]],
    cells: list[str] | None = None,
    chains: dict | None = None,
    progress=None,
    cache: Path | None = TIMING_CACHE,
) -> tuple[dict[str, dict[int, float]], dict[str, Any]]:
    """The same for bins on several hg38 chromosomes, read from local copies of the 11 bigWigs.

    Values already computed are kept per chromosome under `cache` and not read again; only the
    missing bins are lifted and summarised.
    """
    cells = cells or list(REPLI_CELLS)
    t0 = time.time()
    out: dict[str, dict[int, float]] = {}
    missing: dict[str, list[int]] = {}
    known: dict[str, dict[int, float | None]] = {}
    for chrom, bins in wanted.items():
        stored: dict[int, float | None] = {}
        if cache is not None and (cache / f"{chrom}.json.gz").exists():
            with gzip.open(cache / f"{chrom}.json.gz", "rt") as fh:
                stored = {int(k): v for k, v in json.load(fh).items()}
        known[chrom] = stored
        need = [b for b in set(bins) if b not in stored]
        if need:
            missing[chrom] = sorted(need)
    points: dict[tuple[str, int], tuple[str, int]] = {}
    for chrom, bins in missing.items():
        chain = (chains or {}).get(chrom) or load_chain(chrom)
        starts = [c[0] for c in chain]
        for b in bins:
            hit = lift_point(chain, starts, b * BACKGROUND_BIN + BACKGROUND_BIN // 2)
            if hit:
                points[(chrom, b)] = hit
    by_q: dict[str, list[int]] = {}
    for qc, qp in points.values():
        by_q.setdefault(qc, []).append(qp)
    sums: dict[tuple[str, int], float] = {}
    seen: dict[tuple[str, int], int] = {}
    if missing:
        for cell in cells:
            bw = BigWig(repli_local(cell))
            try:
                values: dict[tuple[str, int], float] = {}
                for qc, pts in by_q.items():
                    if qc not in bw.chroms:
                        continue
                    uniq = sorted(set(pts))
                    for st in bw.summarise(qc, [(p, p + 1) for p in uniq], threshold=0.0):
                        if st.bases:
                            values[(qc, st.start)] = st.mean
                for key, hit in points.items():
                    if hit in values:
                        sums[key] = sums.get(key, 0.0) + values[hit]
                        seen[key] = seen.get(key, 0) + 1
            finally:
                bw.close()
            if progress:
                progress(cell)
    for chrom, bins in wanted.items():
        stored = known[chrom]
        for b in missing.get(chrom, []):
            key = (chrom, b)
            stored[b] = sums[key] / seen[key] if seen.get(key) == len(cells) else None
        if cache is not None and chrom in missing:
            cache.mkdir(parents=True, exist_ok=True)
            with gzip.open(cache / f"{chrom}.json.gz", "wt") as fh:
                json.dump(stored, fh)
        out[chrom] = {b: stored[b] for b in bins if stored.get(b) is not None}
    cost = {
        "bins_computed": sum(len(v) for v in missing.values()),
        "bins_cached": sum(len(v) for v in wanted.values()) - sum(len(v) for v in missing.values()),
        "seconds": round(time.time() - t0, 1),
    }
    return out, cost


def rt_edges(rt: dict[int, float], strata: int = RT_STRATA) -> list[float]:
    vals = sorted(rt.values())
    return [vals[int(len(vals) * k / strata)] for k in range(1, strata)] if vals else []


def rt_stratum(value: float | None, edges: list[float]) -> int | None:
    """0 late ... strata-1 early (the signal is high where replication is early)."""
    return None if value is None else bisect_right(edges, value)


# ================================================================================================
# Fast slicing of a chromosome's events
# ================================================================================================
class EventIndex:
    """A chromosome's events, sliced by interval without rescanning the store."""

    def __init__(self, evs: list[Event]):
        self.points = sorted((e for e in evs if e.kind in ("snv", "insertion")), key=lambda e: e.start)
        self.spans = sorted((e for e in evs if e.kind in ("deletion", "structural")), key=lambda e: e.start)
        self.p_starts = [e.start for e in self.points]
        self.s_starts = [e.start for e in self.spans]
        self.max_span = max((e.end - e.start for e in self.spans), default=0)

    def over(self, intervals: list[tuple[int, int]]) -> list[Event]:
        out: list[Event] = []
        seen: set[int] = set()
        for s, e in sorted(intervals):
            out.extend(self.points[bisect_left(self.p_starts, s) : bisect_left(self.p_starts, e)])
            i = bisect_left(self.s_starts, s - self.max_span)
            for ev in self.spans[i : bisect_left(self.s_starts, e)]:
                if ev.end > s and id(ev) not in seen:
                    seen.add(id(ev))
                    out.append(ev)
        out.sort(key=lambda ev: (ev.start, ev.kind))
        return out


def clip_touched(evs: list[Event], intervals: list[tuple[int, int]]) -> tuple[int, int]:
    """Bases of the intervals inside a non-SNV event, and inside a recurring one."""
    touched: set[int] = set()
    recurring: set[int] = set()
    for ev in evs:
        if ev.kind == "snv":
            continue
        rec = ev.minor >= MIN_RECURRING
        for s, e in intervals:
            lo, hi = max(s, ev.start), min(e, ev.end)
            if lo < hi:
                touched.update(range(lo, hi))
                if rec:
                    recurring.update(range(lo, hi))
    return len(touched), len(recurring)


# ================================================================================================
# Structural and population tracks over a set of units
# ================================================================================================
SV_MAX_LENGTH = 1_000_000  # gnomAD's complex calls spanning whole arms are not a unit's structure


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def track_rows(key: str, chrom: str, intervals: list[tuple[int, int]]) -> tuple[list[dict], dict]:
    bb = open_bigbed(key)
    t0 = time.time()
    r0, b0 = bb.src.requests, bb.src.bytes_fetched
    rows = bb.query(chrom, merge_intervals(intervals))
    cost = {
        "requests": bb.src.requests - r0,
        "mb_fetched": round((bb.src.bytes_fetched - b0) / 1e6, 2),
        "seconds": round(time.time() - t0, 1),
    }
    return rows, cost


SV_COMMON_AF = 0.01  # a structural variant carried by 1% of alleles or more is common
HPRC_ARR_ASSEMBLIES = 90  # the arrangement tracks count genomes of the release-1 set


def sv_row_carriers(key: str, row: dict) -> tuple[int | None, int | None, float | None]:
    """(carriers, sample size, frequency) as the track states them."""
    if key.startswith("hprc_") and key != "hprc_v21_sv":
        label = str(row.get("label") or row.get("score") or "")
        try:
            n = int(label.split(":", 1)[0])
        except ValueError:
            return None, HPRC_ARR_ASSEMBLIES, None
        return n, HPRC_ARR_ASSEMBLIES, round(n / HPRC_ARR_ASSEMBLIES, 4)
    if key == "hprc_v21_sv":
        return int(row["AC"]), int(row["alleleNumber"]), float(row["alleleFreq"])
    if key == "gnomad_sv":
        try:
            return int(row["ac"]), int(row["an"]), float(row["af"])
        except (KeyError, ValueError):
            return None, None, None
    return None, None, None


class Overlaps:
    """Rows of one track over a chromosome, asked by interval; a zero-length item (an insertion
    point) is read as the base after it."""

    def __init__(self, rows: list[dict]):
        rows = [
            dict(r, chromEnd=r["chromStart"] + 1) if r["chromEnd"] <= r["chromStart"] else r for r in rows
        ]
        self.rows = sorted(rows, key=lambda r: r["chromStart"])
        self.starts = [r["chromStart"] for r in self.rows]
        self.max_len = max((r["chromEnd"] - r["chromStart"] for r in self.rows), default=0)

    def over(self, s: int, e: int) -> list[dict]:
        i = bisect_left(self.starts, s - self.max_len)
        return [r for r in self.rows[i : bisect_left(self.starts, e)] if r["chromEnd"] > s]


def structure_over(tracks: dict[str, Overlaps], intervals: list[tuple[int, int]]) -> dict[str, Any]:
    """Per structural track: items touching the unit, bases covered by any item and by common ones
    (two or more genomes for the HPRC arrangements, 1% of the track's alleles for the others), and the
    most carriers of any item."""
    bases = sum(e - s for s, e in intervals)
    out: dict[str, Any] = {}
    for key, ov in tracks.items():
        covered: set[int] = set()
        common: set[int] = set()
        n = 0
        top: tuple | None = None
        for s, e in intervals:
            for r in ov.over(s, e):
                if key == "gnomad_sv" and (
                    r.get("FILTER") != "PASS" or r["chromEnd"] - r["chromStart"] > SV_MAX_LENGTH
                ):
                    continue
                n += 1
                span = range(max(s, r["chromStart"]), min(e, r["chromEnd"]))
                covered.update(span)
                c = sv_row_carriers(key, r)
                arrangement = key.startswith("hprc_") and key != "hprc_v21_sv"
                if (arrangement and (c[0] or 0) >= MIN_RECURRING) or (
                    not arrangement and c[2] is not None and c[2] >= SV_COMMON_AF
                ):
                    common.update(span)
                if c[0] is not None and (top is None or c[0] > top[0]):
                    top = c
        if n:
            out[key] = {
                "items": n,
                "covered_share": round(len(covered) / bases, 4) if bases else None,
                "common_covered_share": round(len(common) / bases, 4) if bases else None,
                "most_carriers": top[0] if top else None,
                "sample": top[1] if top else None,
            }
    return out


def gnomad_for_events(rows: Overlaps, evs: list[Event], n: int) -> list[dict[str, Any]]:
    """Each recurring event beside the gnomAD genomes site it corresponds to, the two counts kept apart."""
    out = []
    for ev in evs:
        panel = Counter(t for t in ev.tokens if t is not None)
        m = sum(panel.values())
        entry: dict[str, Any] = {
            "kind": ev.kind,
            "pos": ev.start + 1,
            "length": ev.end - ev.start if ev.kind != "insertion" else None,
            "panel": {
                "haplotypes": m,
                "minor": ev.minor,
                "minor_share": round(ev.minor / m, 3) if m else None,
            },
        }
        match = None
        for r in rows.over(ev.start - 1, ev.end + 1):
            ref, alt = r.get("ref", ""), r.get("alt", "")
            if ev.kind == "snv" and r["chromStart"] == ev.start and len(ref) == len(alt) == 1:
                alts = {t for t in ev.tokens if isinstance(t, str) and t in "ACGT"}
                if alt in alts:
                    match = r
                    break
            elif (
                ev.kind == "deletion"
                and r["chromStart"] == ev.start - 1
                and len(ref) - len(alt) == ev.end - ev.start
            ) or (ev.kind == "insertion" and r["chromStart"] == ev.start and len(alt) > len(ref)):
                match = r
                break
        if match is not None:
            with contextlib.suppress(KeyError, ValueError):
                entry["gnomad"] = {
                    "rsid": match.get("rsId") or None,
                    "af": float(match["AF"]),
                    "an": int(match["AN"]),
                    "ref": match.get("ref"),
                    "alt": match.get("alt"),
                    "filter": match.get("FILTER"),
                }
        out.append(entry)
    _ = n
    return out


# ================================================================================================
# The chromosome: blocks, controls, units, catalogues
# ================================================================================================
UNIT_CLASSES = ("fixed", "storage", "hypervariable", "unplaced")
BLOCK_CLASSES = ("core", "variable", "polymorphic", "lineage_restricted", "unplaced")
TIERS = ("structural", "fossil", "regulatory", "constrained_unknown", "neutral")
MATCHED_WINDOWS = 5  # control windows per block, same length, GC within MATCH_GC and the same timing stratum
MATCH_GC = 0.02
MATCH_TRIES = 400
SEED = 20260913


@functools.lru_cache(maxsize=32)
def coding_genes(chrom: str) -> list[tuple[str, list[tuple[int, int]]]]:
    """Canonical CDS intervals per protein-coding gene from GENCODE (0-based, half-open)."""
    from genomeos.attribution.variation import merge
    from genomeos.genome.annotation import Annotation, default_gencode

    gff = default_gencode({chrom})
    if gff is None:
        return []
    ann = Annotation.from_gff3(gff, {chrom})
    out = []
    for g in ann.protein_coding():
        if g.locus.chrom != chrom:
            continue
        ts = list(g.transcripts.values())
        canon = [t for t in ts if "Ensembl_canonical" in t.tags] or sorted(
            ts, key=lambda t: -(t.locus.end - t.locus.start)
        )
        if canon and canon[0].cds:
            out.append((g.symbol, merge([(c.start, c.end) for c, _ in canon[0].cds])))
    return out


def unit_stats(panel: Panel, index: EventIndex, piece: list[tuple[int, int]]) -> dict[str, Any]:
    d = domain(panel, piece, index.over(piece))
    d.pop("_values")
    d.pop("_events")
    return d


def rt_over(rt: dict[int, float], intervals: list[tuple[int, int]]) -> float | None:
    vals = [
        rt[b]
        for s, e in intervals
        for b in range(s // BACKGROUND_BIN, (e - 1) // BACKGROUND_BIN + 1)
        if b in rt
    ]
    return sum(vals) / len(vals) if vals else None


def class_shares(rows: list[dict], key: str = "class", classes: tuple = UNIT_CLASSES) -> dict[str, Any]:
    c = Counter(r[key] for r in rows)
    n = sum(c.values())
    return {k: {"n": c.get(k, 0), "share": round(c.get(k, 0) / n, 4) if n else None} for k in classes}


def matched_expectation(rows: list[dict], reference: dict[Any, Counter], key) -> dict[str, float | None]:
    """The class shares the reference units would show with these units' strata (placed units only)."""
    exp: Counter = Counter()
    m = 0
    for r in rows:
        if r["class"] == "unplaced":
            continue
        ref = reference.get(key(r))
        if not ref:
            continue
        tot = sum(v for k, v in ref.items() if k != "unplaced")
        if not tot:
            continue
        for k in ("fixed", "storage", "hypervariable"):
            exp[k] += ref.get(k, 0) / tot
        m += 1
    return {k: round(exp[k] / m, 4) if m else None for k in ("fixed", "storage", "hypervariable")}


def placed_shares(rows: list[dict]) -> dict[str, float | None]:
    placed = [r for r in rows if r["class"] != "unplaced"]
    c = Counter(r["class"] for r in placed)
    return {
        k: round(c.get(k, 0) / len(placed), 4) if placed else None
        for k in ("fixed", "storage", "hypervariable")
    }


def spearman(x: list[float], y: list[float]) -> float | None:
    if len(x) < 3:
        return None

    def ranks(v: list[float]) -> list[float]:
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2
            i = j + 1
        return r

    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return round(num / den, 4) if den else None


def block_row(
    panel: Panel,
    index: EventIndex,
    ivs: list[tuple[int, int]],
    bg_gc: Background,
    bg_rt: Background,
    rt: dict[int, float],
    edges: list[float],
) -> dict[str, Any]:
    evs = index.over(ivs)
    m = measure(panel, ivs, evs)
    obs = m["recurring_events"]
    gc_only = block_class(m, obs, bg_gc.expected(panel, ivs))
    with_rt = block_class(m, obs, bg_rt.expected(panel, ivs))
    timing = rt_over(rt, ivs)
    return {
        "measure": m,
        "gc": round(panel.gc_fraction(ivs) or 0, 3) if panel.gc_fraction(ivs) is not None else None,
        "replication_timing": round(timing, 2) if timing is not None else None,
        "rt_stratum": rt_stratum(timing, edges),
        "class_gc_matched": gc_only,
        "class": with_rt,
    }


def matched_windows(
    panel: Panel,
    index: EventIndex,
    target: dict[str, Any],
    length: int,
    forbidden: list[tuple[int, int]],
    bg_gc: Background,
    bg_rt: Background,
    rt: dict[int, float],
    edges: list[float],
    rng,
    k: int = MATCHED_WINDOWS,
) -> list[dict[str, Any]]:
    """Windows of the block's length elsewhere on the chromosome: GC within MATCH_GC, same timing stratum."""
    out: list[dict[str, Any]] = []
    fb = merge_intervals(forbidden)
    fb_starts = [s for s, _ in fb]
    fb_ends = [e for _, e in fb]
    lo, hi = panel.block_start[0], panel.block_end[-1] - length
    if hi <= lo:
        return out
    tries = 0
    while len(out) < k and tries < MATCH_TRIES:
        tries += 1
        s = rng.randrange(lo, hi)
        e = s + length
        i = bisect_right(fb_ends, s)
        if i < len(fb_starts) and fb_starts[i] < e:
            continue
        if panel.aligned_bases(s, e) < 0.9 * length:
            continue
        gc = panel.gc_fraction([(s, e)])
        if gc is None or target["gc"] is None or abs(gc - target["gc"]) > MATCH_GC:
            continue
        st = rt_stratum(rt_over(rt, [(s, e)]), edges)
        if st != target["rt_stratum"]:
            continue
        row = block_row(panel, index, [(s, e)], bg_gc, bg_rt, rt, edges)
        out.append({"start": s, "end": e, "class": row["class"]["class"], "ratio": row["class"]["ratio"]})
    return out


def build(
    chrom: str = "chr21",
    results_dir: Path | None = None,
    cache: Path = CACHE,
    width: int = UNIT_WIDTH,
    progress=None,
    tracks: bool = True,
) -> dict[str, Any]:
    """Read one chromosome's store against its budget, its coding exons and its human axis."""
    import random

    from genomeos.attribution.variation import GNOCCHI_THRESHOLD, gnocchi_over
    from genomeos.results import RESULTS_DIR, load_result

    results_dir = results_dir or RESULTS_DIR
    say = progress or (lambda *a: None)
    t0 = time.time()
    panel = Panel(cache / chrom)
    budget = load_result(f"budget_{chrom}", results_dir) or {}
    variation = load_result(f"variation_{chrom}", results_dir) or {}
    duplication = load_result(f"duplication_{chrom}", results_dir) or {}
    length = budget.get("chromosome_length") or panel.block_end[-1]
    genes = coding_genes(chrom)
    cds_all = merge_intervals([iv for _, ivs in genes for iv in ivs])
    say("events")
    evs = events(panel, [(0, length)])
    index = EventIndex(evs)
    cost: dict[str, Any] = {"maf": {k: panel.meta.get(k) for k in ("mb_fetched", "requests", "seconds")}}

    # replication timing per kilobase, cached beside the store
    rt_path = cache / chrom / "replication_timing.json"
    if rt_path.exists():
        rt = {int(k): v for k, v in json.loads(rt_path.read_text()).items()}
        cost["replication_timing"] = "cached"
    else:
        say("replication timing")
        bins = sorted(
            {
                b
                for i in range(len(panel.block_start))
                for b in range(
                    panel.block_start[i] // BACKGROUND_BIN, (panel.block_end[i] - 1) // BACKGROUND_BIN + 1
                )
            }
        )
        rt, cost["replication_timing"] = replication_timing(chrom, bins)
        rt_path.write_text(json.dumps(rt))
    edges = rt_edges(rt)
    bg_gc = build_background(panel, cds_all, evs)
    bg_rt = build_background(panel, cds_all, evs, rt, edges)

    # -- blocks -------------------------------------------------------------------------------
    say("blocks")
    var_blocks = {(b["start"], b["end"]): b for b in variation.get("blocks", [])}
    dup_blocks = {(b["start"], b["end"]): b for b in duplication.get("blocks", [])}
    blocks = budget.get("blocks", [])
    all_block_ivs = [(b["start"], b["end"]) for b in blocks]
    track_data: dict[str, Overlaps] = {}
    if tracks:
        wanted = all_block_ivs + cds_all
        for key in (
            "hprc_dup",
            "hprc_inv",
            "hprc_del",
            "hprc_ins",
            "hprc_double",
            "hprc_v21_sv",
            "gnomad_sv",
        ):
            rows, c = track_rows(key, chrom, wanted)
            track_data[key] = Overlaps(rows)
            cost[key] = c
    block_rows: list[dict[str, Any]] = []
    for b in blocks:
        ivs = [(b["start"], b["end"])]
        row = block_row(panel, index, ivs, bg_gc, bg_rt, rt, edges)
        v = var_blocks.get((b["start"], b["end"]), {})
        d = dup_blocks.get((b["start"], b["end"]), {})
        row.update(
            {
                "start": b["start"],
                "end": b["end"],
                "length": b["length"],
                "sequence_class": b["class"],
                "tier": b["guess"]["tier"],
                "mammal_fraction": (b.get("phylop") or {}).get("fraction_above"),
                "gnocchi_case": (v.get("case") or {}).get("case") or "unmeasured",
                "human_fraction": (v.get("gnocchi") or {}).get("fraction_above"),
                "duplicated_fraction": d.get("duplicated_fraction"),
                "structure": structure_over(track_data, ivs) if track_data else {},
            }
        )
        block_rows.append(row)

    # -- the coding control, gene by gene ------------------------------------------------------
    say("coding control")
    gene_rows = []
    for symbol, ivs in genes:
        row = block_row(panel, index, ivs, bg_gc, bg_rt, rt, edges)
        row.update({"gene": symbol, "cds_bases": sum(e - s for s, e in ivs)})
        gene_rows.append(row)
    pooled = {}
    for name, ivs in [("cds", cds_all)] + [
        (t, [(b["start"], b["end"]) for b in blocks if b["guess"]["tier"] == t]) for t in TIERS
    ]:
        if not ivs:
            continue
        m = measure(panel, ivs, index.over(ivs))
        pooled[name] = {
            "measure": m,
            "ratio_gc": round(m["recurring_events"] / bg_gc.expected(panel, ivs), 3),
            "ratio_gc_rt": round(m["recurring_events"] / bg_rt.expected(panel, ivs), 3),
        }

    def tally(rows: list[dict], key: str) -> dict[str, int]:
        c = Counter(r[key]["class"] for r in rows)
        return {k: c.get(k, 0) for k in BLOCK_CLASSES}

    def core_share(counts: dict[str, int]) -> float | None:
        callable_ = counts["core"] + counts["variable"]
        return round(counts["core"] / callable_, 3) if callable_ else None

    by_tier = {}
    for t in TIERS:
        rows = [r for r in block_rows if r["tier"] == t]
        by_tier[t] = {
            "blocks": len(rows),
            "bp": sum(r["length"] for r in rows),
            "classes": tally(rows, "class"),
            "classes_gc_only": tally(rows, "class_gc_matched"),
        }
        by_tier[t]["core_share_of_callable"] = core_share(by_tier[t]["classes"])
        by_tier[t]["core_share_of_callable_gc_only"] = core_share(by_tier[t]["classes_gc_only"])
    genes_tally = {
        "genes": len(gene_rows),
        "classes": tally(gene_rows, "class"),
        "classes_gc_only": tally(gene_rows, "class_gc_matched"),
    }
    genes_tally["core_share_of_callable"] = core_share(genes_tally["classes"])
    genes_tally["core_share_of_callable_gc_only"] = core_share(genes_tally["classes_gc_only"])

    # -- matched control windows for every callable block --------------------------------------
    say("matched windows")
    rng = random.Random(SEED)
    forbidden = cds_all + [(b["start"], b["end"]) for b in blocks if b["guess"]["tier"] == "structural"]
    matched_by_tier: dict[str, Counter] = {t: Counter() for t in TIERS}
    matched_cds: Counter = Counter()
    for r in block_rows:
        if r["class"]["class"] == "unplaced" or r["tier"] == "structural" or r["length"] > 200_000:
            continue
        for w in matched_windows(panel, index, r, r["length"], forbidden, bg_gc, bg_rt, rt, edges, rng, k=2):
            matched_by_tier[r["tier"]][w["class"]] += 1
    for r in gene_rows:
        if r["class"]["class"] == "unplaced":
            continue
        span = sum(
            e - s for s, e in [(iv[0], iv[1]) for iv in genes[[g for g, _ in genes].index(r["gene"])][1]]
        )
        for w in matched_windows(panel, index, r, span, forbidden, bg_gc, bg_rt, rt, edges, rng, k=2):
            matched_cds[w["class"]] += 1
    matched = {
        "cds_genes": {k: matched_cds.get(k, 0) for k in BLOCK_CLASSES},
        **{t: {k: v.get(k, 0) for k in BLOCK_CLASSES} for t, v in matched_by_tier.items()},
    }
    for v in matched.values():
        v["core_share_of_callable"] = core_share(v)

    # -- units ----------------------------------------------------------------------------------
    say("units")

    def unit_rows_for(
        group: str, pieces: list[list[tuple[int, int]]], block_id: int | None = None
    ) -> list[dict]:
        out = []
        for piece in pieces:
            d = unit_stats(panel, index, piece)
            timing = rt_over(rt, piece)
            d.update(
                {
                    "group": group,
                    "block": block_id,
                    "start": piece[0][0],
                    "end": piece[-1][1],
                    "gc_stratum": gc_stratum(panel.gc_fraction(piece)),
                    "rt_stratum": rt_stratum(timing, edges),
                }
            )
            out.append(d)
        return out

    units: list[dict] = []
    for symbol, ivs in genes:
        units += unit_rows_for("cds", tile(ivs, width))
        _ = symbol
    for i, b in enumerate(blocks):
        if b["class"] == "gap":
            continue
        units += unit_rows_for(b["guess"]["tier"], tile([(b["start"], b["end"])], width), i)
    unknown_ivs = merge_intervals([(b["start"], b["end"]) for b in blocks])
    background_pieces = []
    cds_starts = [s for s, _ in cds_all]
    cds_ends = [e for _, e in cds_all]
    grid_lo = (panel.block_start[0] // width) * width
    for s in range(grid_lo, panel.block_end[-1], width):
        e = s + width
        i = bisect_right(cds_ends, s)
        if i < len(cds_starts) and cds_starts[i] < e:
            continue
        if panel.aligned_bases(s, e) < width * UNIT_MIN_ALIGNED:
            continue
        background_pieces.append([(s, e)])
    background = unit_rows_for("background", background_pieces)
    ref_gc: dict[Any, Counter] = {}
    ref_gcrt: dict[Any, Counter] = {}
    for u in background:
        ref_gc.setdefault(u["gc_stratum"], Counter())[u["class"]] += 1
        ref_gcrt.setdefault((u["gc_stratum"], u["rt_stratum"]), Counter())[u["class"]] += 1
    unit_groups = {}
    for group in ("cds", *TIERS, "background"):
        rows = background if group == "background" else [u for u in units if u["group"] == group]
        if not rows:
            continue
        unit_groups[group] = {
            "units": len(rows),
            "bp": sum(u["bases"] for u in rows),
            "classes": class_shares(rows),
            "placed_shares": placed_shares(rows),
            "expected_gc_matched": matched_expectation(rows, ref_gc, lambda r: r["gc_stratum"]),
            "expected_gc_rt_matched": matched_expectation(
                rows, ref_gcrt, lambda r: (r["gc_stratum"], r["rt_stratum"])
            ),
        }
    # timing of the unit classes, chromosome-wide
    rt_by_class = {}
    for k in ("fixed", "storage", "hypervariable", "unplaced"):
        rows = [u for u in background if u["class"] == k]
        c = Counter(u["rt_stratum"] for u in rows)
        n = sum(v for s, v in c.items() if s is not None)
        rt_by_class[k] = {
            ("late", "middle", "early")[s]: round(c.get(s, 0) / n, 3) if n else None for s in range(RT_STRATA)
        }
        rt_by_class[k]["units"] = len(rows)

    result: dict[str, Any] = {
        "chrom": chrom,
        "assemblies": panel.n,
        "panel": {
            "reference": REFERENCE,
            "track": MAF_TRACK,
            "source": panel.meta.get("source"),
            "sites": len(panel.site_pos),
            "blocks": len(panel.block_start),
        },
        "evidence": SOURCE_NOTES | {"replication_timing": RT_EVIDENCE},
        "thresholds": thresholds(width),
        "background": {
            "recurring_per_kb_overall": round(1000 * bg_gc.overall, 3),
            "recurring_per_kb_by_gc_stratum": {
                str(k): round(1000 * v, 3)
                for k, v in sorted(bg_gc.density.items(), key=lambda kv: (kv[0] is None, kv[0]))
            },
            "recurring_per_kb_by_gc_and_timing": {
                f"gc{k[0]}_rt{k[1]}": round(1000 * v, 3)
                for k, v in sorted(bg_rt.density.items(), key=lambda kv: str(kv[0]))
            },
            "rt_edges": [round(x, 2) for x in edges],
        },
        "pooled": pooled,
        "blocks_by_tier": by_tier,
        "genes": genes_tally,
        "matched_windows": matched,
        "units": unit_groups,
        "timing_by_unit_class": rt_by_class,
    }
    result["_block_rows"] = block_rows
    result["_gene_rows"] = gene_rows
    result["_units"] = units
    result["_background_units"] = background
    result["_index"] = index
    result["_panel"] = panel
    result["_rt"] = rt
    result["_edges"] = edges
    result["_bg"] = (bg_gc, bg_rt)
    result["_genes"] = genes
    result["_unknown_ivs"] = unknown_ivs
    result["cost"] = cost
    result["seconds_core"] = round(time.time() - t0, 1)
    _ = (GNOCCHI_THRESHOLD, gnocchi_over)
    return result


def thresholds(width: int = UNIT_WIDTH) -> dict[str, Any]:
    return {
        "min_recurring": MIN_RECURRING,
        "unit_width": width,
        "unit_min_aligned": UNIT_MIN_ALIGNED,
        "unit_min_informative": UNIT_MIN_INFORMATIVE,
        "fixed_min_share": FIXED_MIN_SHARE,
        "storage_min_cover": STORAGE_MIN_COVER,
        "storage_max_values": STORAGE_MAX_VALUES,
        "core_max_ratio": CORE_MAX_RATIO,
        "core_max_p": CORE_MAX_P,
        "min_expected": MIN_EXPECTED,
        "present_min": PRESENT_MIN,
        "lineage_max": LINEAGE_MAX,
        "structural_max": STRUCTURAL_MAX,
        "aligned_min": ALIGNED_MIN,
        "missing_max": MISSING_MAX,
        "gc_strata": list(GC_STRATA),
        "rt_strata": RT_STRATA,
        "matched_gc": MATCH_GC,
    }


# ================================================================================================
# The storage catalogue: value domains, population frequencies, the executor shortlist
# ================================================================================================
def storage_catalogue(result: dict[str, Any], chrom: str, examples: int = 12) -> dict[str, Any]:
    """Every storage unit of the UNKNOWN space with its value domain, frequencies and read-out evidence."""
    panel: Panel = result["_panel"]
    index: EventIndex = result["_index"]
    units = [u for u in result["_units"] if u["group"] in TIERS and u["class"] == "storage"]
    hyper = [u for u in result["_units"] if u["group"] in TIERS and u["class"] == "hypervariable"]
    ivs = [(u["start"], u["end"]) for u in units + hyper]
    cost: dict[str, Any] = {}
    tables: dict[str, Overlaps] = {}
    for key in ("gnomad_snv", "trexplorer", "gtex_dapg", "mpravardb"):
        rows, cost[key] = track_rows(key, chrom, ivs)
        if key == "gtex_dapg":  # an eQTL row spans variant to gene; keep the variant's base
            rows = [
                dict(r, chromStart=r_pos, chromEnd=r_pos + 1)
                for r in rows
                for r_pos in [_eqtl_pos(r)]
                if r_pos is not None
            ]
        tables[key] = Overlaps(rows)
    catalogue = []
    kinds: Counter = Counter()
    values_hist: Counter = Counter()
    matched = total_snv_indel = 0
    panel_vs_gnomad: list[tuple[float, float]] = []
    for u in units:
        piece = [(u["start"], u["end"])]
        d = domain(panel, piece, index.over(piece))
        rec: list[Event] = d["_events"]
        kind = "+".join(sorted({ev.kind for ev in rec}))
        kinds[kind] += 1
        values_hist[min(d["recurring_values"], 9)] += 1
        freq = gnomad_for_events(tables["gnomad_snv"], rec, panel.n)
        for f in freq:
            if f["kind"] in ("snv", "deletion", "insertion"):
                total_snv_indel += 1
                if "gnomad" in f:
                    matched += 1
                    af = f["gnomad"]["af"]
                    if f["panel"]["minor_share"] is not None:
                        panel_vs_gnomad.append((f["panel"]["minor_share"], min(af, 1 - af)))
        trs = [
            {
                "motif": r.get("referenceMotif"),
                "motif_size": int(r.get("motifSize") or 0),
                "ref_copies": r.get("numRepeats"),
                "tenk_hist": r.get("tenKAlleleHist") or None,
                "tenk_alleles": int(r.get("tenKNumAlleles") or 0),
                "hprc_hist": r.get("hprcAlleleHist") or None,
            }
            for r in tables["trexplorer"].over(u["start"], u["end"])
            if int(r.get("motifSize") or 0) >= 2
        ]
        eqtl = sorted(
            {(r.get("eqtlName"), r.get("geneName")) for r in tables["gtex_dapg"].over(u["start"], u["end"])}
        )
        mpra = [
            {
                "rsid": r.get("rsid"),
                "cell": r.get("cellLine"),
                "log2fc": r.get("log2FC"),
                "fdr": r.get("fdr"),
                "study": (r.get("mpraStudy") or "")[:80],
            }
            for r in tables["mpravardb"].over(u["start"], u["end"])
        ]
        entry = {
            "start": u["start"],
            "end": u["end"],
            "tier": u["group"],
            "block": u["block"],
            "rt_stratum": u["rt_stratum"],
            "domain": {
                k: d[k]
                for k in (
                    "informative",
                    "recurring_events",
                    "values",
                    "recurring_values",
                    "top_share",
                    "recurring_cover",
                    "entropy_bits",
                    "effective_values",
                )
            },
            "kinds": kind,
            "values": describe_values(d, 6),
            "length_alleles": length_domain(panel, piece, rec)
            if any(ev.kind != "snv" for ev in rec)
            else None,
            "frequencies": freq,
            "tandem_repeats": trs,
            "eqtl": [{"variant": a, "gene": b} for a, b in eqtl],
            "mpra": mpra,
        }
        catalogue.append(entry)
    hyper_by_length: Counter = Counter()
    for u in hyper:
        piece = [(u["start"], u["end"])]
        rec = [ev for ev in index.over(piece) if ev.minor >= MIN_RECURRING]
        if any(ev.kind != "snv" for ev in rec):
            hyper_by_length[length_class(length_domain(panel, piece, rec), panel.n)["class"]] += 1
        else:
            hyper_by_length["substitutions only"] += 1
    with_eqtl = [c for c in catalogue if c["eqtl"]]
    with_mpra = [c for c in catalogue if c["mpra"]]
    tr_units = [c for c in catalogue if c["tandem_repeats"]]
    hyper_tr = sum(1 for u in hyper if tables["trexplorer"].over(u["start"], u["end"]))
    summary = {
        "storage_units": len(catalogue),
        "storage_bp": sum(c["end"] - c["start"] for c in catalogue),
        "by_tier": dict(Counter(c["tier"] for c in catalogue)),
        "recurring_values_histogram": {
            ("9+" if k == 9 else str(k)): v for k, v in sorted(values_hist.items())
        },
        "effective_values_median": _median([c["domain"]["effective_values"] for c in catalogue]),
        "top_share_median": _median([c["domain"]["top_share"] for c in catalogue]),
        "event_kinds": dict(kinds.most_common(8)),
        "gnomad_matched_events": matched,
        "panel_events_snv_indel": total_snv_indel,
        "panel_minor_share_vs_gnomad_maf_spearman": spearman(
            [a for a, _ in panel_vs_gnomad], [b for _, b in panel_vs_gnomad]
        ),
        "with_tandem_repeat": len(tr_units),
        "with_gtex_eqtl": len(with_eqtl),
        "with_mpra_allele_pair": len(with_mpra),
        "with_either": len([c for c in catalogue if c["eqtl"] or c["mpra"]]),
        "hypervariable_units": len(hyper),
        "hypervariable_with_tandem_repeat": hyper_tr,
        "hypervariable_read_by_length": dict(hyper_by_length),
        "sample_sizes": {
            "panel": f"{panel.n} haplotype assemblies (hg38 is the coordinate system and not counted)",
            "gnomad": "gnomAD v4.1.1 genomes, AN per site (up to 152,430 alleles)",
            "trexplorer": "TenK10K allele counts per locus as stated (typically 3,850), HPRC where given",
        },
    }

    def pick(rows: list[dict], n: int) -> list[dict]:
        return sorted(
            rows,
            key=lambda c: (
                -c["domain"]["informative"],
                c["domain"]["recurring_values"],
                -c["domain"]["entropy_bits"],
            ),
        )[:n]

    showcase = {
        "executor_shortlist": pick([c for c in catalogue if c["eqtl"] or c["mpra"]], examples),
        "constrained_unknown": pick([c for c in catalogue if c["tier"] == "constrained_unknown"], examples),
        "two_value_slots": pick(
            [c for c in catalogue if c["domain"]["recurring_values"] == 2], examples // 2
        ),
        "tandem_repeat_lengths": pick(tr_units, examples // 2),
    }
    for rows in showcase.values():
        for c in rows:
            c["frequencies"] = [f for f in c["frequencies"] if f["panel"]["minor"] >= MIN_RECURRING][:8]
    return {"summary": summary, "showcase": showcase, "cost": cost, "_catalogue": catalogue}


def _eqtl_pos(row: dict) -> int | None:
    m = re.match(r"(\w+):(\d+)", row.get("eqtlPos") or "")
    return int(m.group(2)) - 1 if m else None


def _median(vals: list[float | None]) -> float | None:
    v = sorted(x for x in vals if x is not None)
    if not v:
        return None
    k = len(v) // 2
    return v[k] if len(v) % 2 else round((v[k - 1] + v[k]) / 2, 3)


# ================================================================================================
# The human axis beside Gnocchi, kilobase by kilobase and block by block
# ================================================================================================
DEPLETED_P = 0.05  # a kilobase is depleted of recurring variation when its Poisson lower tail is this small


def against_gnocchi(result: dict[str, Any], chrom: str, results_dir: Path | None = None) -> dict[str, Any]:
    from genomeos.attribution.variation import GNOCCHI_THRESHOLD, gnocchi_over
    from genomeos.results import RESULTS_DIR

    results_dir = results_dir or RESULTS_DIR
    panel: Panel = result["_panel"]
    index: EventIndex = result["_index"]
    bg_rt: Background = result["_bg"][1]
    edges = result["_edges"]
    rt = result["_rt"]
    unknown = result["_unknown_ivs"]
    u_starts = [s for s, _ in unknown]
    u_ends = [e for _, e in unknown]
    dups = _superdups(chrom, results_dir)
    d_starts = [s for s, _ in dups]
    d_ends = [e for _, e in dups]
    bins = sorted(b for b, row in bg_rt.bins.items() if row[0] >= 0.9 * BACKGROUND_BIN and row[2] is not None)
    stats, cost = gnocchi_over(chrom, [(b * BACKGROUND_BIN, (b + 1) * BACKGROUND_BIN) for b in bins])
    rows = []
    for b, st in zip(bins, stats, strict=True):
        s, e = b * BACKGROUND_BIN, (b + 1) * BACKGROUND_BIN
        row = bg_rt.bins[b]
        exp = row[0] * bg_rt.rate(b)
        obs = row[1]
        below, _ = poisson_tails(int(obs), exp)
        evs = index.over([(s, e)])
        rec = [ev for ev in evs if ev.minor >= MIN_RECURRING]
        m = measure(panel, [(s, e)], evs)
        i = bisect_right(u_ends, s)
        j = bisect_right(d_ends, s)
        rows.append(
            {
                "bin": b,
                "z": st.mean if st.bases else None,
                "ratio": obs / exp if exp else None,
                "depleted": below <= DEPLETED_P,
                "indel_share": (sum(ev.kind != "snv" for ev in rec) / len(rec)) if rec else None,
                "touched_recurring": m["touched_recurring_share"],
                "missing": m["missing_share"],
                "rt": rt_stratum(rt.get(b), edges),
                "gc": row[2],
                "unknown": i < len(u_starts) and u_starts[i] < e,
                "duplicated": j < len(d_starts) and d_starts[j] < e,
            }
        )
    scored = [r for r in rows if r["z"] is not None and r["ratio"] is not None]
    groups = {
        "both": [r for r in scored if r["z"] >= GNOCCHI_THRESHOLD and r["depleted"]],
        "gnocchi_only": [r for r in scored if r["z"] >= GNOCCHI_THRESHOLD and not r["depleted"]],
        "panel_only": [r for r in scored if r["z"] < GNOCCHI_THRESHOLD and r["depleted"]],
        "neither": [r for r in scored if r["z"] < GNOCCHI_THRESHOLD and not r["depleted"]],
        "gnocchi_unscored_panel_depleted": [r for r in rows if r["z"] is None and r["depleted"]],
        "gnocchi_unscored_other": [r for r in rows if r["z"] is None and not r["depleted"]],
    }

    def describe(rs: list[dict]) -> dict[str, Any]:
        if not rs:
            return {"kilobases": 0}
        c = Counter(r["rt"] for r in rs)
        n = sum(v for k, v in c.items() if k is not None)
        return {
            "kilobases": len(rs),
            "mean_ratio": round(sum(r["ratio"] for r in rs if r["ratio"] is not None) / len(rs), 3),
            "mean_z": round(
                sum(r["z"] for r in rs if r["z"] is not None) / max(1, sum(r["z"] is not None for r in rs)), 2
            )
            if any(r["z"] is not None for r in rs)
            else None,
            "indel_share_of_recurring": _median([r["indel_share"] for r in rs]),
            "mean_touched_recurring": round(sum(r["touched_recurring"] or 0 for r in rs) / len(rs), 4),
            "mean_missing": round(sum(r["missing"] or 0 for r in rs) / len(rs), 4),
            "duplicated_share": round(sum(r["duplicated"] for r in rs) / len(rs), 3),
            "unknown_share": round(sum(r["unknown"] for r in rs) / len(rs), 3),
            "mean_gc": round(sum(r["gc"] for r in rs) / len(rs), 3),
            "timing": {
                ("late", "middle", "early")[k]: round(c.get(k, 0) / n, 3) if n else None
                for k in range(RT_STRATA)
            },
        }

    per_kb = {k: describe(v) for k, v in groups.items()}
    per_kb["spearman_z_vs_ratio"] = spearman([r["z"] for r in scored], [r["ratio"] for r in scored])
    per_kb["kilobases_read"] = len(rows)
    per_kb["definition"] = (
        "a kilobase outside canonical CDS with 90% or more aligned; "
        f"Gnocchi constrained at mean Z >= {GNOCCHI_THRESHOLD}; "
        f"panel depleted at a Poisson lower tail <= {DEPLETED_P} against the GC and timing matched rate"
    )
    # blocks: the panel's class against the case the two older axes made
    table: dict[str, Counter] = {}
    for r in result["_block_rows"]:
        table.setdefault(r["class"]["class"], Counter())[r["gnocchi_case"]] += 1
    cases = ("syntax", "relaxed", "recent", "tolerant", "unmeasured")
    blocks = {cls: {c: table.get(cls, Counter()).get(c, 0) for c in cases} for cls in BLOCK_CLASSES}
    human = Counter()
    for r in result["_block_rows"]:
        hf = r["human_fraction"]
        cls = r["class"]["class"]
        if hf is None or cls not in ("core", "variable"):
            continue
        human[(cls, hf >= 0.25)] += 1
    return {
        "per_kilobase": per_kb,
        "blocks_class_by_case": blocks,
        "blocks_core_or_variable_by_gnocchi": {
            "core_gnocchi_constrained": human[("core", True)],
            "core_gnocchi_free": human[("core", False)],
            "variable_gnocchi_constrained": human[("variable", True)],
            "variable_gnocchi_free": human[("variable", False)],
        },
        "cost": cost,
        "_rows": rows,
    }


def depletion_null(result: dict[str, Any], p: float = DEPLETED_P) -> dict[str, Any]:
    """How many kilobases a Poisson would call depleted by chance, and how overdispersed the counts are.

    Recurring events share genealogies (a haplotype carries many of them at once), so counts per
    kilobase vary far more than a Poisson allows; the dispersion says by how much, and why the
    matched control windows, not the Poisson tail, are the test for a block.
    """
    bg: Background = result["_bg"][1]
    n = observed = 0
    expected = 0.0
    z: list[float] = []
    for b, row in bg.bins.items():
        if row[0] < 0.9 * BACKGROUND_BIN or row[2] is None:
            continue
        e = row[0] * bg.rate(b)
        if e <= 0:
            continue
        n += 1
        cdf = 0.0
        term = math.exp(-e)
        k = 0
        while cdf + term <= p:
            cdf += term
            k += 1
            term *= e / k
        expected += cdf
        observed += poisson_tails(int(row[1]), e)[0] <= p
        z.append((row[1] - e) / math.sqrt(e))
    mean = sum(z) / len(z) if z else 0.0
    return {
        "kilobases": n,
        "depleted_observed": observed,
        "depleted_expected_under_poisson": round(expected),
        "pearson_dispersion": round(sum((x - mean) ** 2 for x in z) / len(z), 2) if z else None,
    }


def _superdups(chrom: str, results_dir: Path) -> list[tuple[int, int]]:
    p = results_dir / f"superdups_{chrom}.bed.gz"
    if not p.exists():
        return []
    ivs = []
    with gzip.open(p, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2 and f[0].isdigit():
                ivs.append((int(f[0]), int(f[1])))
    return merge_intervals(ivs)


# ================================================================================================
# Sensitivity: how the classes move with the thresholds
# ================================================================================================
def sensitivity(result: dict[str, Any], widths: tuple[int, ...] = (100, 500, 1000)) -> dict[str, Any]:
    panel: Panel = result["_panel"]
    index: EventIndex = result["_index"]
    units = result["_units"]
    background = result["_background_units"]
    out: dict[str, Any] = {"unit_thresholds": {}, "unit_width": {}, "block_thresholds": {}}
    groups = {
        "cds": [u for u in units if u["group"] == "cds"],
        "neutral": [u for u in units if u["group"] == "neutral"],
        "background": background,
    }
    for fixed_min in (0.9, 0.95, 0.99):
        for max_values in (4, 8, 16):
            key = f"fixed>={fixed_min},storage<={max_values}"
            out["unit_thresholds"][key] = {
                g: _shares_with(rows, panel.n, fixed_min, STORAGE_MIN_COVER, max_values)
                for g, rows in groups.items()
            }
    for w in widths:
        row = {}
        for g, pieces in (
            ("cds", [p for _, ivs in result["_genes"] for p in tile(ivs, w)]),
            (
                "neutral",
                [
                    p
                    for r in result["_block_rows"]
                    if r["tier"] == "neutral"
                    for p in tile([(r["start"], r["end"])], w)
                ],
            ),
        ):
            rows = [unit_stats(panel, index, p) for p in pieces]
            row[g] = placed_shares([dict(r, **{"class": r["class"]}) for r in rows])
        out["unit_width"][str(w)] = row
    bg_rt = result["_bg"][1]
    for ratio in (0.3, 0.5, 0.7):
        for present in (0.9, 0.95, 0.99):
            counts = {}
            for name, rows in (
                ("cds_genes", result["_gene_rows"]),
                ("neutral", [r for r in result["_block_rows"] if r["tier"] == "neutral"]),
                (
                    "constrained_unknown",
                    [r for r in result["_block_rows"] if r["tier"] == "constrained_unknown"],
                ),
                ("regulatory", [r for r in result["_block_rows"] if r["tier"] == "regulatory"]),
            ):
                c = Counter(_reclass(r, ratio, present) for r in rows)
                callable_ = c["core"] + c["variable"]
                counts[name] = {
                    "core": c["core"],
                    "variable": c["variable"],
                    "core_share_of_callable": round(c["core"] / callable_, 3) if callable_ else None,
                }
            out["block_thresholds"][f"ratio<={ratio},presence>={present}"] = counts
    _ = bg_rt
    return out


def _shares_with(
    rows: list[dict], n: int, fixed_min: float, cover_min: float, max_values: int
) -> dict[str, float | None]:
    cls = [unit_class(r, n, fixed_min, cover_min, max_values) for r in rows]
    placed = [c for c in cls if c != "unplaced"]
    k = Counter(placed)
    return {
        c: round(k.get(c, 0) / len(placed), 4) if placed else None
        for c in ("fixed", "storage", "hypervariable")
    }


def _reclass(r: dict, core_ratio: float, present_min: float) -> str:
    c = r["class"]
    m = r["measure"]
    if c["class"] in ("unplaced", "lineage_restricted"):
        return c["class"]
    if (m["presence"] or 0) < present_min or (m["touched_recurring_share"] or 0) >= STRUCTURAL_MAX:
        return "polymorphic"
    if c["expected"] < MIN_EXPECTED:
        return "unplaced"
    if c["ratio"] is not None and c["ratio"] <= core_ratio and c["p_depleted"] <= CORE_MAX_P:
        return "core"
    return "variable"


# ================================================================================================
# Calibration: slots whose values are known, read with the same instrument
# ================================================================================================
CALIBRATION_FLANK = 20_000  # the regional stores hold this much either side, the local background
CALIBRATION = (
    {
        "name": "HERC2 intron 86, the OCA2 enhancer (rs12913832)",
        "chrom": "chr15",
        "kind": "variant",
        "pos0": 28120471,
        "rsid": "rs12913832",
        "expect": "two recurring values (A brown, G blue)",
    },
    {
        "name": "MCM6 intron 13, the LCT enhancer (rs4988235)",
        "chrom": "chr2",
        "kind": "variant",
        "pos0": 135851075,
        "rsid": "rs4988235",
        "expect": "two recurring values (lactase persistence)",
    },
    {
        "name": "ABO coding exons",
        "chrom": "chr9",
        "kind": "cds",
        "gene": "ABO",
        "expect": "a handful of values (A, B, O alleles; hg38 carries O, rs8176719)",
    },
    {
        "name": "HLA-A coding exons",
        "chrom": "chr6",
        "kind": "cds",
        "gene": "HLA-A",
        "expect": "many values: the most polymorphic human gene",
    },
    {
        "name": "INS promoter VNTR",
        "chrom": "chr11",
        "kind": "tandem_repeat",
        "gene": "INS",
        "motif_size": 14,
        "expect": "length classes (class I short, class III long)",
    },
    {
        "name": "DRD4 exon 3 VNTR",
        "chrom": "chr11",
        "kind": "tandem_repeat",
        "gene": "DRD4",
        "motif_size": 48,
        "expect": "repeat counts 2, 4, 7 recurring",
    },
    {
        "name": "SLC6A3 3' UTR VNTR",
        "chrom": "chr5",
        "kind": "tandem_repeat",
        "gene": "SLC6A3",
        "motif_size": 40,
        "expect": "repeat counts 9 and 10 recurring",
    },
)
VARIANT_HALF_WIDTH = 100


def _gene_span(chrom: str, symbol: str) -> tuple[int, int, list[tuple[int, int]]] | None:
    for sym, ivs in coding_genes(chrom):
        if sym == symbol:
            return ivs[0][0], ivs[-1][1], ivs
    return None


def calibration_units(loci=CALIBRATION) -> list[dict[str, Any]]:
    """Resolve each calibration locus to hg38 intervals, verifying looked-up positions against the sources."""
    out = []
    genes_cache: dict[str, list] = {}
    for loc in loci:
        chrom = loc["chrom"]
        row = dict(loc)
        if loc["kind"] == "variant":
            rows, _ = track_rows("gnomad_snv", chrom, [(loc["pos0"], loc["pos0"] + 1)])
            hit = [r for r in rows if r.get("rsId") == loc["rsid"] and r["chromStart"] == loc["pos0"]]
            row["verified"] = bool(hit)
            row["gnomad"] = (
                {
                    "ref": hit[0]["ref"],
                    "alt": hit[0]["alt"],
                    "af": float(hit[0]["AF"]),
                    "an": int(hit[0]["AN"]),
                }
                if hit
                else None
            )
            row["intervals"] = [(loc["pos0"] - VARIANT_HALF_WIDTH, loc["pos0"] + VARIANT_HALF_WIDTH)]
        elif loc["kind"] == "cds":
            genes = genes_cache.setdefault(chrom, coding_genes(chrom))
            ivs = next((iv for sym, iv in genes if sym == loc["gene"]), None)
            row["verified"] = ivs is not None
            row["intervals"] = ivs or []
        else:
            genes = genes_cache.setdefault(chrom, coding_genes(chrom))
            ivs = next((iv for sym, iv in genes if sym == loc["gene"]), None)
            span = (ivs[0][0] - 5_000, ivs[-1][1] + 5_000) if ivs else None
            best = None
            if span:
                rows, _ = track_rows("trexplorer", chrom, [span])
                cands = [r for r in rows if int(r.get("motifSize") or 0) == loc["motif_size"]]
                best = max(cands, key=lambda r: r["chromEnd"] - r["chromStart"], default=None)
            row["verified"] = best is not None
            if best:
                row["intervals"] = [(best["chromStart"] - 50, best["chromEnd"] + 50)]
                row["trexplorer"] = {
                    "locus": [best["chromStart"], best["chromEnd"]],
                    "motif": best.get("referenceMotif"),
                    "ref_copies": best.get("numRepeats"),
                    "tenk_hist": best.get("tenKAlleleHist") or None,
                    "tenk_alleles": int(best.get("tenKNumAlleles") or 0),
                    "hprc_hist": best.get("hprcAlleleHist") or None,
                    "hprc_alleles": int(best.get("hprcNumAlleles") or 0),
                }
            else:
                row["intervals"] = []
        out.append(row)
    return out


def read_locus(
    panel: Panel,
    intervals: list[tuple[int, int]],
    flank: int = CALIBRATION_FLANK,
    rt_value: float | None = None,
    edges=None,
) -> dict[str, Any]:
    """A unit in a regional store: its domain, its values, and its rate against the flanks around it."""
    evs = events(panel, intervals)
    d = domain(panel, intervals, evs)
    lo, hi = intervals[0][0], intervals[-1][1]
    flanks = [(max(0, lo - flank), lo), (hi, hi + flank)]
    fev = events(panel, flanks)
    fm = measure(panel, flanks, fev)
    um = measure(panel, intervals, evs)
    density = fm["recurring_events"] / fm["aligned_bases"] if fm["aligned_bases"] else 0.0
    exp = density * um["aligned_bases"]
    tiles = [domain(panel, p) for p in tile(flanks, UNIT_WIDTH)]
    return {
        "bases": d["bases"],
        "class": d["class"],
        "domain": {
            k: d[k]
            for k in (
                "informative",
                "recurring_events",
                "values",
                "recurring_values",
                "top_share",
                "recurring_cover",
                "entropy_bits",
                "effective_values",
            )
        },
        "values": describe_values(d, 10),
        "alternatives": alternatives(d)[:12],
        "length_alleles": length_domain(panel, intervals, d["_events"]),
        "length_class": length_class(length_domain(panel, intervals, d["_events"]), panel.n),
        "class_by_max_values": {str(k): unit_class(d, panel.n, max_values=k) for k in (4, 8, 16, 32)},
        "measure": um,
        "against_flanks": {
            "recurring_events": um["recurring_events"],
            "expected_from_flanks": round(exp, 1),
            "ratio": round(um["recurring_events"] / exp, 3) if exp else None,
            "flank_units": placed_shares(tiles),
        },
        "replication_timing": round(rt_value, 2) if rt_value is not None else None,
        "rt_stratum_on_chr21_scale": rt_stratum(rt_value, edges) if edges else None,
        "_events": d["_events"],
    }


def calibrate(cache: Path = CACHE, edges: list[float] | None = None) -> list[dict[str, Any]]:
    loci = calibration_units()
    wanted: dict[str, list[int]] = {}
    for loc in loci:
        for s, e in loc["intervals"]:
            wanted.setdefault(loc["chrom"], []).extend(
                range(s // BACKGROUND_BIN, (e - 1) // BACKGROUND_BIN + 1)
            )
    timing, _ = replication_timing_many(wanted)
    out = []
    for loc in loci:
        row = dict(loc)
        store = cache / f"{loc['chrom']}_regions"
        if not loc["intervals"] or not (store / "meta.json").exists():
            row["read"] = None
            out.append(row)
            continue
        panel = Panel(store)
        rt = timing.get(loc["chrom"], {})
        read = read_locus(panel, loc["intervals"], rt_value=rt_over(rt, loc["intervals"]), edges=edges)
        frows, _ = track_rows("gnomad_snv", loc["chrom"], loc["intervals"])
        read["frequencies"] = gnomad_for_events(Overlaps(frows), read.pop("_events"), panel.n)
        row["read"] = read
        out.append(row)
    return out


# ================================================================================================
# The 69 candidates, each against its own flanks
# ================================================================================================
def read_candidates(
    results_dir: Path | None = None, cache: Path = CACHE, edges: list[float] | None = None
) -> dict[str, Any]:
    from genomeos.results import RESULTS_DIR, load_result

    results_dir = results_dir or RESULTS_DIR
    cands = (load_result("syntax_candidates_genome_wide", results_dir) or {}).get("candidates", [])
    rows = []
    by_chrom: dict[str, list[dict]] = {}
    for c in cands:
        by_chrom.setdefault(c["chrom"], []).append(c)
    timing, _ = replication_timing_many(
        {
            chrom: sorted(
                {
                    b
                    for c in cs
                    for b in range(c["start"] // BACKGROUND_BIN, (c["end"] - 1) // BACKGROUND_BIN + 1)
                }
            )
            for chrom, cs in by_chrom.items()
        }
    )
    for chrom, cs in sorted(by_chrom.items()):
        store = cache / f"{chrom}_regions"
        if not (store / "meta.json").exists():
            continue
        panel = Panel(store)
        rt = timing.get(chrom, {})
        for c in cs:
            ivs = [(c["start"], c["end"])]
            read = read_locus(panel, ivs, rt_value=rt_over(rt, ivs), edges=edges)
            read.pop("_events")
            m = read["measure"]
            obs = read["against_flanks"]["recurring_events"]
            exp = read["against_flanks"]["expected_from_flanks"]
            cls = block_class(m, obs, exp)
            tiles = [domain(panel, p) for p in tile(ivs, UNIT_WIDTH)]
            rows.append(
                {
                    "chrom": chrom,
                    "start": c["start"],
                    "end": c["end"],
                    "length": c["length"],
                    "reading": (c.get("reading") or {}).get("class"),
                    "mammal_fraction": c.get("mammal_fraction"),
                    "human_fraction": c.get("human_fraction"),
                    "class": cls,
                    "presence": m["presence"],
                    "missing_share": m["missing_share"],
                    "pi": m["pi"],
                    "touched_recurring_share": m["touched_recurring_share"],
                    "units": dict(Counter(t["class"] for t in tiles)),
                    "replication_timing": read["replication_timing"],
                    "rt_stratum_on_chr21_scale": read["rt_stratum_on_chr21_scale"],
                }
            )
    c = Counter(r["class"]["class"] for r in rows)
    by_reading: dict[str, Counter] = {}
    for r in rows:
        by_reading.setdefault(r["reading"] or "none", Counter())[r["class"]["class"]] += 1
    ratios = [r["class"]["ratio"] for r in rows if r["class"]["ratio"] is not None]
    return {
        "read": len(rows),
        "of": len(cands),
        "classes": {k: c.get(k, 0) for k in BLOCK_CLASSES},
        "by_reading": {k: dict(v) for k, v in by_reading.items()},
        "median_ratio_to_flanks": _median(ratios),
        "pooled_ratio_to_flanks": round(
            sum(r["class"]["recurring_events"] for r in rows)
            / max(1e-9, sum(r["class"]["expected"] for r in rows)),
            3,
        )
        if rows
        else None,
        "rows": rows,
    }


def local_candidates(result: dict[str, Any], chrom: str, results_dir: Path | None = None) -> dict[str, Any]:
    """The organiser's real unknown of this chromosome (constrained_unknown blocks that are not
    copies), each read in the whole-chromosome store against its own flanks, as the 69 were."""
    from genomeos.results import RESULTS_DIR, load_result

    results_dir = results_dir or RESULTS_DIR
    panel: Panel = result["_panel"]
    rt = result["_rt"]
    edges = result["_edges"]
    blocks = (load_result(f"organised_{chrom}", results_dir) or {}).get("candidates", [])
    rows = []
    for b in blocks:
        ivs = [(b["start"], b["end"])]
        read = read_locus(panel, ivs, rt_value=rt_over(rt, ivs), edges=edges)
        read.pop("_events")
        cls = block_class(
            read["measure"],
            read["against_flanks"]["recurring_events"],
            read["against_flanks"]["expected_from_flanks"],
        )
        rows.append(
            {
                "start": b["start"],
                "end": b["end"],
                "length": b["length"],
                "case": b.get("case"),
                "mammal_fraction": b.get("mammal_fraction"),
                "human_fraction": b.get("human_fraction"),
                "class": cls,
                "presence": read["measure"]["presence"],
                "missing_share": read["measure"]["missing_share"],
                "domain_class": read["class"],
                "replication_timing": read["replication_timing"],
            }
        )
    obs = sum(r["class"]["recurring_events"] for r in rows)
    exp = sum(r["class"]["expected"] for r in rows)
    return {
        "source": f"organised_{chrom} candidates (constrained_unknown, not copies)",
        "read": len(rows),
        "classes": dict(Counter(r["class"]["class"] for r in rows)),
        "pooled_ratio_to_flanks": round(obs / exp, 3) if exp else None,
        "median_ratio_to_flanks": _median([r["class"]["ratio"] for r in rows]),
        "rows": rows,
    }


# ================================================================================================
# What the genome would cost
# ================================================================================================
MAIN_CHROMOSOMES = tuple(f"chr{i}" for i in [*range(1, 23), "X", "Y"])


def genome_cost(result: dict[str, Any], storage: dict[str, Any] | None = None) -> dict[str, Any]:
    sizes = {}
    for c in MAIN_CHROMOSOMES:
        try:
            sizes[c] = remote_size(MAF_URL.format(chrom=c))
        except OSError:
            sizes[c] = None
    total = sum(v for v in sizes.values() if v)
    maf = result["cost"]["maf"]
    rate = maf["mb_fetched"] / maf["seconds"] if maf.get("seconds") else None
    store_mb = sum(p.stat().st_size for p in (CACHE / result["chrom"]).glob("*")) / 1e6
    chr21_mb = sizes.get(result["chrom"]) or 1
    scale = total / chr21_mb
    out = {
        "maf_gb_by_chromosome": {k: round(v / 1e9, 2) if v else None for k, v in sizes.items()},
        "maf_gb_total": round(total / 1e9, 1),
        "measured_mb_per_second": round(rate, 1) if rate else None,
        "maf_stream_hours": round(total / 1e6 / rate / 3600, 1) if rate else None,
        "store_mb_chr21": round(store_mb, 1),
        "store_gb_genome_estimate": round(store_mb * scale / 1e3, 2),
        "memory_note": (
            "a store is held in memory with its events; chr21 used about 1 GB, "
            "chr1 and chr2 would need 5 to 6 GB"
        ),
        "scale_factor_from_chr21": round(scale, 1),
    }
    if storage:
        mb = storage["cost"].get("gnomad_snv", {}).get("mb_fetched")
        if mb:
            out["gnomad_mb_chr21_storage_units"] = mb
            out["gnomad_gb_genome_estimate"] = round(mb * scale / 1e3, 1)
    return out


# ================================================================================================
# The committed summary
# ================================================================================================
def compact_block(r: dict[str, Any], units: Counter | None = None) -> dict[str, Any]:
    m = r["measure"]
    c = r["class"]
    g = r["class_gc_matched"]
    return {
        **{k: r[k] for k in ("start", "end", "length", "tier", "sequence_class") if k in r},
        **({"gene": r["gene"], "cds_bases": r["cds_bases"]} if "gene" in r else {}),
        "class": c["class"],
        "confidence": c["confidence"],
        "reason": c.get("reason"),
        "recurring_events": c["recurring_events"],
        "expected": c["expected"],
        "ratio": c["ratio"],
        "p_depleted": c["p_depleted"],
        "class_gc_only": g["class"],
        "ratio_gc_only": g["ratio"],
        "presence": m["presence"],
        "deleted_share": m["deleted_share"],
        "replaced_share": m["replaced_share"],
        "missing_share": m["missing_share"],
        "aligned_share": m["aligned_share"],
        "assemblies_aligned_90": m["assemblies_aligned_90"],
        "identity": m["identity"],
        "pi": m["pi"],
        "touched_share": m["touched_share"],
        "touched_recurring_share": m["touched_recurring_share"],
        "gc": r["gc"],
        "replication_timing": r["replication_timing"],
        "rt_stratum": r["rt_stratum"],
        **({"units": dict(units)} if units is not None else {}),
        **{
            k: r[k]
            for k in ("gnocchi_case", "human_fraction", "mammal_fraction", "duplicated_fraction")
            if k in r
        },
        **({"structure": r["structure"]} if r.get("structure") else {}),
    }


STRUCTURE_TRACKS = ("hprc_dup", "hprc_inv", "hprc_del", "hprc_ins", "hprc_double", "hprc_v21_sv", "gnomad_sv")


def structure_by_tier(blocks: list[dict[str, Any]]) -> dict[str, Any]:
    """Bases of each tier under any item and under common items of every structural track."""
    out: dict[str, Any] = {}
    for t in TIERS:
        rows = [b for b in blocks if b["tier"] == t]
        bp = sum(b["length"] for b in rows) or 1
        out[t] = {
            key: {
                "any": round(
                    sum(
                        (b.get("structure", {}).get(key, {}).get("covered_share") or 0) * b["length"]
                        for b in rows
                    )
                    / bp,
                    4,
                ),
                "common": round(
                    sum(
                        (b.get("structure", {}).get(key, {}).get("common_covered_share") or 0) * b["length"]
                        for b in rows
                    )
                    / bp,
                    4,
                ),
            }
            for key in STRUCTURE_TRACKS
        }
    return out


def refresh_structure(chrom: str, results_dir: Path | None = None) -> dict[str, Any]:
    """Re-read the structural tracks over a saved result's blocks and rewrite those fields only."""
    from genomeos.results import RESULTS_DIR, load_result, save_result

    results_dir = results_dir or RESULTS_DIR
    res = load_result(f"human_panel_{chrom}", results_dir)
    ivs = [(b["start"], b["end"]) for b in res["blocks"]]
    tracks: dict[str, Overlaps] = {}
    cost = {}
    for key in STRUCTURE_TRACKS:
        rows, cost[key] = track_rows(key, chrom, ivs)
        tracks[key] = Overlaps(rows)
    for b in res["blocks"]:
        b["structure"] = structure_over(tracks, [(b["start"], b["end"])])
    res["structure_by_tier"] = structure_by_tier(res["blocks"])
    res["structure_definition"] = (
        f"any item, and recurring or common items: HPRC arrangements carried by {MIN_RECURRING} or more of "
        f"{HPRC_ARR_ASSEMBLIES} genomes, v2.1 SVs and gnomAD SVs (PASS, under {SV_MAX_LENGTH:,} bp) at AF >= "
        f"{SV_COMMON_AF}; the block classes use the alignment's own events, these tracks are annotation"
    )
    res["cost"].update({f"refresh_{k}": v for k, v in cost.items()})
    payload = {k: v for k, v in res.items() if k not in ("result", "date")}
    save_result(f"human_panel_{chrom}", payload, results_dir)
    return res["structure_by_tier"]


def summarise(
    result: dict[str, Any],
    catalogue: dict[str, Any],
    gnocchi: dict[str, Any],
    sens: dict[str, Any],
    calibration: list[dict[str, Any]],
    candidates: dict[str, Any],
    cost: dict[str, Any],
) -> dict[str, Any]:
    per_block_units: dict[int, Counter] = {}
    for u in result["_units"]:
        if u["block"] is not None:
            per_block_units.setdefault(u["block"], Counter())[u["class"]] += 1
    blocks = [
        compact_block(r, per_block_units.get(i, Counter())) for i, r in enumerate(result["_block_rows"])
    ]
    genes = [compact_block(r) for r in result["_gene_rows"]]
    units = result["units"]
    control = control_verdict(result)
    catalogues = three_catalogues(result)
    out = {k: v for k, v in result.items() if not k.startswith("_") and k not in ("cost", "seconds_core")}
    out.update(
        {
            "control": control,
            "catalogues": catalogues,
            "storage": catalogue["summary"],
            "storage_showcase": catalogue["showcase"],
            "against_gnocchi": {k: v for k, v in gnocchi.items() if not k.startswith("_")},
            "sensitivity": sens,
            "calibration": calibration,
            "candidates": candidates,
            "genome_wide": cost,
            "cost": {**result["cost"], "storage_tracks": catalogue["cost"], "gnocchi": gnocchi.get("cost")},
            "structure_by_tier": structure_by_tier(blocks),
            "blocks": blocks,
            "coding_genes": genes,
        }
    )
    _ = units
    return out


def control_verdict(result: dict[str, Any]) -> dict[str, Any]:
    """Does the instrument see coding exons as held and the neutral tier as free? Read before the biology."""
    u = result["units"]
    genes = result["genes"]
    neutral = result["blocks_by_tier"]["neutral"]
    pooled = result["pooled"]
    cds_fixed = u["cds"]["placed_shares"]["fixed"]
    neu_fixed = u["neutral"]["placed_shares"]["fixed"]
    checks = {
        "cds_units_fixed_above_matched": cds_fixed - u["cds"]["expected_gc_rt_matched"]["fixed"] >= 0.15,
        "neutral_units_fixed_at_matched": abs(neu_fixed - u["neutral"]["expected_gc_rt_matched"]["fixed"])
        <= 0.05,
        "cds_pooled_ratio_below_0.6": pooled["cds"]["ratio_gc_rt"] <= 0.6,
        "neutral_pooled_ratio_0.8_to_1.25": 0.8 <= pooled["neutral"]["ratio_gc_rt"] <= 1.25,
        "cds_gene_core_share_exceeds_neutral_by_0.3": (genes["core_share_of_callable"] or 0)
        - (neutral["core_share_of_callable"] or 0)
        >= 0.3,
        "cds_gene_core_share_exceeds_matched_windows_by_0.3": (genes["core_share_of_callable"] or 0)
        - (result["matched_windows"]["cds_genes"]["core_share_of_callable"] or 0)
        >= 0.3,
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "overwhelmingly_core": (genes["core_share_of_callable"] or 0) >= 0.8,
        "numbers": {
            "cds_units_fixed": cds_fixed,
            "cds_units_fixed_expected_gc_rt": u["cds"]["expected_gc_rt_matched"]["fixed"],
            "neutral_units_fixed": neu_fixed,
            "neutral_units_fixed_expected_gc_rt": u["neutral"]["expected_gc_rt_matched"]["fixed"],
            "cds_pooled_ratio_gc_rt": pooled["cds"]["ratio_gc_rt"],
            "neutral_pooled_ratio_gc_rt": pooled["neutral"]["ratio_gc_rt"],
            "cds_gene_core_share_of_callable": genes["core_share_of_callable"],
            "cds_gene_core_share_gc_only": genes["core_share_of_callable_gc_only"],
            "neutral_block_core_share_of_callable": neutral["core_share_of_callable"],
            "matched_windows_core_share_cds_like": result["matched_windows"]["cds_genes"][
                "core_share_of_callable"
            ],
            "matched_windows_core_share_neutral_like": result["matched_windows"]["neutral"][
                "core_share_of_callable"
            ],
        },
        "rule": (
            "the checks are stated in code (control_verdict) and were written after the first chr21 read "
            "of the "
            "coding and neutral numbers, so they describe the ordering seen rather than a prediction; "
            "'overwhelmingly core' "
            "(80% of callable genes) is reported beside them, not required"
        ),
    }


def three_catalogues(result: dict[str, Any]) -> dict[str, Any]:
    """Albert's three lists for the UNKNOWN space and for coding exons: fixed, storage, cannot place."""
    out = {}
    for name, groups in (("unknown_space", TIERS), ("coding_exons", ("cds",))):
        rows = [u for u in result["_units"] if u["group"] in groups]
        c: Counter = Counter()
        bp: Counter = Counter()
        for u in rows:
            c[u["class"]] += 1
            bp[u["class"]] += u["bases"]
        out[name] = {
            "fixed": {"units": c["fixed"], "bp": bp["fixed"]},
            "storage": {"units": c["storage"], "bp": bp["storage"]},
            "cannot_place": {
                "units": c["hypervariable"] + c["unplaced"],
                "bp": bp["hypervariable"] + bp["unplaced"],
                "hypervariable": {"units": c["hypervariable"], "bp": bp["hypervariable"]},
                "not_aligned_or_missing": {"units": c["unplaced"], "bp": bp["unplaced"]},
            },
        }
    return out


def save_catalogue(catalogue: dict[str, Any], chrom: str, cache: Path = CACHE) -> Path:
    p = cache / chrom / "storage_catalogue.json.gz"
    with gzip.open(p, "wt") as fh:
        json.dump(catalogue["_catalogue"], fh)
    return p
