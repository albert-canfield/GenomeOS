# SPDX-License-Identifier: AGPL-3.0-or-later
"""Measured contact frequency from a 4D Nucleome .hic matrix, read by HTTP range, never downloaded.

The enhancer-to-gene models so far used 1/distance as the contact term: a power law standing in for
the measurement. 4DN publishes the matrices themselves, but a released in-situ Hi-C .hic of a deeply
sequenced cell line is 6 to 23 GB, which no laptop should hold for a few thousand cells of the
matrix. The .hic container is built for exactly that: a header, a master index of compressed blocks
and a normalisation index, all addressable, so the bins a question needs can be read with HTTP range
requests and nothing else. This module is that reader, in the standard library only (struct, zlib,
urllib), for version 8 and version 9 files.

What it answers is one question per call: the normalised contact between two positions on one
chromosome at a chosen resolution. Normalisation is the file's own balancing vector, the first of
SCALE, KR, VC_SQRT and VC that the file carries and that defines both bins of the question, and
observed over expected uses the file's expected-value vector for the same resolution, so a contact
1 Mb away is not compared with a contact 5 kb away on the raw scale. Contacts are `experimental`
evidence: the number is a measurement of a cell population, not a model output.

The blocks are fetched, read and dropped; only the answers are cached, under
data/knowledge/hic_contact (a few MB per cell line), with the cell line, the file's accession and
the resolution in the name, so the cache can never mix two matrices.
"""

from __future__ import annotations

import base64
import json
import math
import struct
import urllib.request
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.genome.hic import FOURDN, HOW, _DropAuthOnRedirect, credentials

KNOWLEDGE = Path("data/knowledge/hic_contact")
EVIDENCE = "experimental: 4D Nucleome in-situ Hi-C contact matrix (.hic), read by range request"
RESOLUTION = 5_000  # the bin size the pairs are asked at; 4DN's merged .hic files carry 1 kb upward
NORM_PREFERENCE = ("SCALE", "KR", "VC_SQRT", "VC")  # balancing vectors, best first
UNIT = "BP"
CHUNK = 1 << 22  # 4 MB per range request while streaming the footer
MAX_EXPECTED_DISTANCE_BINS = 10_000  # beyond this the expected vector's last value is reused
BLOCK_CACHE = 8  # decoded blocks kept in memory; a block covers megabases, so a few suffice

# The cell lines whose matrices the CRISPRi comparison uses, and the 4DN file chosen for each:
# both are "merged replicates" in-situ Hi-C (MboI, Aiden lab), the deepest released GRCh38 matrix
# of that biosource on the portal.
MATRICES = {
    "K562": "4DNFITUOMFUQ",
    "GM12878": "4DNFI1UEG1HD",
}


def download_url(accession: str) -> str:
    return f"{FOURDN}/files-processed/{accession}/@@download/{accession}.hic"


class RangeReader:
    """Bytes of a remote file by HTTP range, with 4DN's key and its redirect to S3."""

    def __init__(self, url: str, auth: tuple[str, str] | None = None, timeout: int = 300) -> None:
        self.url = url
        self.auth = auth
        self.timeout = timeout
        self.requests = 0
        self.bytes = 0

    def read(self, offset: int, length: int) -> bytes:
        if length <= 0:
            return b""
        headers = {"Range": f"bytes={offset}-{offset + length - 1}", "User-Agent": "GenomeOS/0.1 (range)"}
        if self.auth:
            token = base64.b64encode(f"{self.auth[0]}:{self.auth[1]}".encode()).decode()
            headers["Authorization"] = f"Basic {token}"
        req = urllib.request.Request(self.url, headers=headers)
        opener = urllib.request.build_opener(_DropAuthOnRedirect())
        with opener.open(req, timeout=self.timeout) as r:  # noqa: S310
            if r.status != 206:
                raise OSError(f"{self.url}: range request answered {r.status}, not 206 Partial Content")
            data = r.read()
        self.requests += 1
        self.bytes += len(data)
        return data


class _Cursor:
    """A little-endian cursor over a remote file, refilling itself by range request."""

    def __init__(self, reader: RangeReader, offset: int, chunk: int = CHUNK) -> None:
        self.reader = reader
        self.offset = offset
        self.chunk = chunk
        self.buf = b""
        self.pos = 0

    def _need(self, n: int) -> None:
        if self.pos + n <= len(self.buf):
            return
        self.buf = self.buf[self.pos :]
        self.offset += self.pos
        self.pos = 0
        while len(self.buf) < n:
            self.buf += self.reader.read(self.offset + len(self.buf), max(self.chunk, n - len(self.buf)))

    def take(self, n: int) -> bytes:
        self._need(n)
        out = self.buf[self.pos : self.pos + n]
        self.pos += n
        return out

    def skip(self, n: int) -> None:
        if self.pos + n <= len(self.buf):
            self.pos += n
        else:
            self.offset += self.pos + n
            self.pos = 0
            self.buf = b""

    def string(self) -> str:
        out = bytearray()
        while True:
            c = self.take(1)
            if c in (b"\x00", b""):
                return out.decode("utf-8", "replace")
            out += c

    def int32(self) -> int:
        return struct.unpack("<i", self.take(4))[0]

    def int64(self) -> int:
        return struct.unpack("<q", self.take(8))[0]

    def float32(self) -> float:
        return struct.unpack("<f", self.take(4))[0]

    def float64(self) -> float:
        return struct.unpack("<d", self.take(8))[0]


def _chrom_key(name: str) -> str:
    """4DN's .hic files name chromosomes without the prefix; the benchmark tables use chr1."""
    n = name[3:] if name.lower().startswith("chr") else name
    return n.upper() if n.lower() in ("x", "y", "m", "mt") else n


@dataclass
class _Resolution:
    binsize: int
    block_bin_count: int
    block_column_count: int
    blocks: dict[int, tuple[int, int]]  # block number -> (position, compressed size)


class HicMatrix:
    """One remote .hic file: its header, the blocks of one chromosome, and normalised contacts.

    The layout read here is the documented .hic container (versions 8 and 9): header, per-chromosome
    matrix records with a block index, zlib-compressed blocks of (bin1, bin2, count), an
    expected-value vector per resolution and a normalisation vector per chromosome and resolution.
    """

    def __init__(self, url: str, auth: tuple[str, str] | None = None, reader: RangeReader | None = None):
        self.reader = reader or RangeReader(url, auth)
        self.version = 0
        self.footer_position = 0
        self.chroms: list[tuple[str, int]] = []
        self.resolutions: list[int] = []
        self.master: dict[str, tuple[int, int]] = {}
        self.norm_index: dict[tuple[str, int, str, int], tuple[int, int]] = {}
        self.expected: dict[tuple[str, int], tuple[list[float], dict[int, float]]] = {}
        self.norm_index_position = 0
        self._matrices: dict[tuple[int, int], _Resolution | None] = {}
        self._norms: dict[tuple[str, int, int], list[float]] = {}
        self._blocks: dict[tuple[int, int, int], dict[tuple[int, int], float]] = {}
        self._block_order: list[tuple[int, int, int]] = []
        self._read_header()

    # ---- header and footer -------------------------------------------------

    def _read_header(self) -> None:
        c = _Cursor(self.reader, 0, chunk=1 << 16)
        magic = c.string()
        if magic != "HIC":
            raise OSError(f"not a .hic file: magic {magic!r}")
        self.version = c.int32()
        if self.version < 8:
            raise OSError(f".hic version {self.version} is older than the documented v8 layout")
        self.footer_position = c.int64()
        c.string()  # genome id
        if self.version >= 9:
            self.norm_index_position = c.int64()
            c.int64()  # length of the normalisation index, unused: it is parsed entry by entry
        for _ in range(c.int32()):
            c.string(), c.string()  # attributes
        for _ in range(c.int32()):
            name = c.string()
            size = c.int64() if self.version >= 9 else c.int32()
            self.chroms.append((_chrom_key(name), size))
        self.resolutions = [c.int32() for _ in range(c.int32())]

    def chrom_index(self, chrom: str) -> int:
        key = _chrom_key(chrom)
        for i, (name, _) in enumerate(self.chroms):
            if name == key:
                return i
        raise LookupError(f"chromosome {chrom!r} is not in the matrix")

    def chrom_length(self, chrom: str) -> int:
        return self.chroms[self.chrom_index(chrom)][1]

    def read_footer(self, binsize: int = RESOLUTION) -> None:
        """Master index, the expected-value vectors of one resolution, and the whole norm index."""
        if self.master:
            return
        c = _Cursor(self.reader, self.footer_position)
        c.int64() if self.version >= 9 else c.int32()  # bytes of the v5 footer, unused
        for _ in range(c.int32()):
            key = c.string()
            self.master[key] = (c.int64(), c.int32())
        self._read_expected(c, binsize, normalised=False)
        self._read_expected(c, binsize, normalised=True)
        if self.version >= 9 and self.norm_index_position:
            c = _Cursor(self.reader, self.norm_index_position)
        for _ in range(c.int32()):
            kind, chrom, unit, size = c.string(), c.int32(), c.string(), c.int32()
            position = c.int64()
            nbytes = c.int64() if self.version >= 9 else c.int32()
            self.norm_index[(kind, chrom, unit, size)] = (position, nbytes)

    def _read_expected(self, c: _Cursor, binsize: int, normalised: bool) -> None:
        for _ in range(c.int32()):
            kind = c.string() if normalised else ""
            unit, size = c.string(), c.int32()
            n = c.int64() if self.version >= 9 else c.int32()
            width = 4 if self.version >= 9 else 8
            keep = unit == UNIT and size == binsize and (not normalised or kind in NORM_PREFERENCE)
            if keep:
                raw = c.take(n * width)
                values = list(struct.unpack(f"<{n}{'f' if self.version >= 9 else 'd'}", raw))
            else:
                c.skip(n * width)
                values = []
            factors: dict[int, float] = {}
            for _ in range(c.int32()):
                chrom = c.int32()
                factor = c.float32() if self.version >= 9 else c.float64()
                factors[chrom] = factor
            if keep:
                self.expected[(kind or "NONE", size)] = (values, factors)

    # ---- normalisation ----------------------------------------------------

    def norm_types(self, chrom_index: int, binsize: int) -> list[str]:
        """The balancing vectors this file carries for one chromosome and resolution, best first."""
        return [k for k in NORM_PREFERENCE if (k, chrom_index, UNIT, binsize) in self.norm_index]

    def norm_type(self, chrom_index: int, binsize: int) -> str | None:
        kinds = self.norm_types(chrom_index, binsize)
        return kinds[0] if kinds else None

    def norm_vector(self, kind: str, chrom_index: int, binsize: int) -> list[float]:
        key = (kind, chrom_index, binsize)
        if key not in self._norms:
            position, nbytes = self.norm_index[(kind, chrom_index, UNIT, binsize)]
            c = _Cursor(self.reader, position, chunk=max(CHUNK, nbytes))
            n = c.int64() if self.version >= 9 else c.int32()
            width = 4 if self.version >= 9 else 8
            raw = c.take(n * width)
            self._norms[key] = list(struct.unpack(f"<{n}{'f' if self.version >= 9 else 'd'}", raw))
        return self._norms[key]

    def expected_at(self, kind: str, binsize: int, chrom_index: int, distance_bins: int) -> float | None:
        entry = self.expected.get((kind, binsize)) or self.expected.get(("NONE", binsize))
        if not entry or not entry[0]:
            return None
        values, factors = entry
        d = min(distance_bins, len(values) - 1, MAX_EXPECTED_DISTANCE_BINS)
        v = values[d]
        factor = factors.get(chrom_index, 1.0) or 1.0
        return v / factor if v > 0 else None

    # ---- blocks -----------------------------------------------------------

    def _resolution(self, chrom_index: int, binsize: int) -> _Resolution | None:
        key = (chrom_index, binsize)
        if key in self._matrices:
            return self._matrices[key]
        self.read_footer(binsize)
        entry = self.master.get(f"{chrom_index}_{chrom_index}")
        if entry is None:
            self._matrices[key] = None
            return None
        c = _Cursor(self.reader, entry[0], chunk=max(CHUNK, entry[1]))
        c.int32(), c.int32()  # the two chromosome indices
        found = None
        for _ in range(c.int32()):
            unit = c.string()
            c.int32()  # resolution index
            for _ in range(4):
                c.float32()  # sum, occupied cells, standard deviation, 95th percentile
            size = c.int32()
            block_bin_count, block_column_count = c.int32(), c.int32()
            blocks = {}
            for _ in range(c.int32()):
                number = c.int32()
                blocks[number] = (c.int64(), c.int32())
            if unit == UNIT and size == binsize:
                found = _Resolution(size, block_bin_count, block_column_count, blocks)
        self._matrices[key] = found
        return found

    def _block_numbers(self, res: _Resolution, bin1: int, bin2: int) -> list[int]:
        """The blocks that can hold one cell, both triangles, for the v8 and the v9 layouts."""
        x, y = min(bin1, bin2), max(bin1, bin2)
        if self.version < 9:
            r, c = y // res.block_bin_count, x // res.block_bin_count
            return [r * res.block_column_count + c, c * res.block_column_count + r]
        # v9 intra blocks are indexed on the rotated (diagonal, depth) grid
        position = (x + y) // 2 // res.block_bin_count
        depth = int(math.log2(1 + abs(y - x) / math.sqrt(2) / res.block_bin_count))
        return [depth * res.block_column_count + position]

    def _read_block(self, chrom_index: int, binsize: int, number: int, position: int, size: int):
        key = (chrom_index, binsize, number)
        cached = self._blocks.get(key)
        if cached is not None:
            return cached
        raw = zlib.decompress(self.reader.read(position, size))
        records = _decode_block(raw, self.version)
        self._blocks[key] = records
        self._block_order.append(key)
        while len(self._block_order) > BLOCK_CACHE:
            self._blocks.pop(self._block_order.pop(0), None)
        return records

    # ---- the question -----------------------------------------------------

    def contact(self, chrom: str, pos1: int, pos2: int, binsize: int = RESOLUTION) -> dict[str, Any] | None:
        """Normalised contact between the bins of two positions, with observed over expected.

        Balancing uses the first vector of NORM_PREFERENCE that has a finite positive weight for both
        bins: matrix balancing leaves single bins undefined (KR does not converge everywhere), and
        dropping those pairs would silently thin the comparison, while mixing vectors inside one pair
        would not be one measurement. Which vector was used is reported with the number.

        Returns None when the chromosome or the resolution is missing, or when no balancing vector
        defines both bins (a bin no normalisation reaches, typically unmappable).
        """
        ci = self.chrom_index(chrom)
        res = self._resolution(ci, binsize)
        if res is None:
            return None
        b1, b2 = pos1 // binsize, pos2 // binsize
        kind, w1, w2 = "", 0.0, 0.0
        for candidate in self.norm_types(ci, binsize):
            weights = self.norm_vector(candidate, ci, binsize)
            if max(b1, b2) >= len(weights):
                continue
            a, b = weights[b1], weights[b2]
            if a > 0 and b > 0 and not math.isnan(a) and not math.isnan(b):
                kind, w1, w2 = candidate, a, b
                break
        if not kind:
            return None
        raw = 0.0
        seen = False
        for number in dict.fromkeys(self._block_numbers(res, b1, b2)):
            entry = res.blocks.get(number)
            if entry is None:
                continue
            records = self._read_block(ci, binsize, number, *entry)
            for cell in ((b1, b2), (b2, b1)):
                v = records.get(cell)
                if v is not None:
                    raw, seen = v, True
                    break
            if seen:
                break
        observed = raw / (w1 * w2)
        expected = self.expected_at(kind, binsize, ci, abs(b2 - b1))
        return {
            "norm": kind,
            "binsize": binsize,
            "same_bin": b1 == b2,
            "raw": raw,
            "observed": observed,
            "expected": expected,
            "oe": (observed / expected) if expected else None,
        }


def _decode_block(raw: bytes, version: int) -> dict[tuple[int, int], float]:
    """One decompressed block of a .hic matrix: {(bin1, bin2): count}, v8 and v9 layouts."""
    out: dict[tuple[int, int], float] = {}
    n_records = struct.unpack_from("<i", raw, 0)[0]
    p = 4
    if version < 7:  # pragma: no cover - refused in the header
        raise OSError("block layout older than v7")
    bin_x_offset, bin_y_offset = struct.unpack_from("<ii", raw, p)
    p += 8
    use_short = raw[p] == 0
    p += 1
    short_x = short_y = True
    if version >= 9:
        short_x, short_y = raw[p] == 0, raw[p + 1] == 0
        p += 2
    kind = raw[p]
    p += 1

    def take(fmt: str, width: int):
        nonlocal p
        v = struct.unpack_from(fmt, raw, p)[0]
        p += width
        return v

    if kind == 1:
        rows = take("<h", 2) if short_y else take("<i", 4)
        for _ in range(rows):
            y = bin_y_offset + (take("<h", 2) if short_y else take("<i", 4))
            cols = take("<h", 2) if short_x else take("<i", 4)
            for _ in range(cols):
                x = bin_x_offset + (take("<h", 2) if short_x else take("<i", 4))
                value = take("<h", 2) if use_short else take("<f", 4)
                out[(x, y)] = float(value)
    elif kind == 2:
        n_points = take("<i", 4)
        width = take("<h", 2)
        for i in range(n_points):
            value = take("<h", 2) if use_short else take("<f", 4)
            if use_short and value == -32768:
                continue
            if not use_short and math.isnan(value):
                continue
            row, col = divmod(i, width)
            out[(bin_x_offset + col, bin_y_offset + row)] = float(value)
    else:
        raise OSError(f"unknown block layout {kind}")
    if n_records and not out and kind == 1:  # pragma: no cover - a layout we misread would be silent
        raise OSError("block decoded to no records")
    return out


class ContactSource:
    """Contacts for several cell lines, cached on disk, each from its own 4DN matrix.

    `contact(cell, chrom, pos1, pos2)` answers from the cache when it can and reads the matrix when
    it cannot. Nothing but the answers is written: the blocks are decompressed in memory and dropped.
    """

    def __init__(
        self,
        matrices: dict[str, str] | None = None,
        knowledge: Path = KNOWLEDGE,
        binsize: int = RESOLUTION,
        auth: tuple[str, str] | None = None,
        readers: dict[str, HicMatrix] | None = None,
    ) -> None:
        self.matrices = dict(matrices or MATRICES)
        self.knowledge = knowledge
        self.binsize = binsize
        self.auth = auth
        self._open: dict[str, HicMatrix | None] = dict(readers or {})
        self._cache: dict[str, dict[str, Any]] = {}
        self._dirty: set[str] = set()
        self.misses = 0

    # cache file per cell line, resolution and accession: two matrices can never mix
    def _path(self, cell: str) -> Path:
        return self.knowledge / f"{cell}_{self.matrices[cell]}_{self.binsize}.json"

    def _load(self, cell: str) -> dict[str, Any]:
        if cell not in self._cache:
            p = self._path(cell)
            self._cache[cell] = json.loads(p.read_text()) if p.exists() else {}
        return self._cache[cell]

    def _matrix(self, cell: str) -> HicMatrix | None:
        if cell not in self._open:
            accession = self.matrices.get(cell)
            if accession is None:
                self._open[cell] = None
            else:
                auth = self.auth or credentials()
                if auth is None:
                    raise PermissionError(f"reading a 4DN matrix needs an account key: {HOW}")
                self._open[cell] = HicMatrix(download_url(accession), auth)
        return self._open[cell]

    def contact(self, cell: str, chrom: str, pos1: int, pos2: int) -> dict[str, Any] | None:
        if cell not in self.matrices:
            return None
        cache = self._load(cell)
        b1, b2 = pos1 // self.binsize, pos2 // self.binsize
        key = f"{_chrom_key(chrom)}:{min(b1, b2)}:{max(b1, b2)}"
        if key in cache:
            return cache[key]
        self.misses += 1
        m = self._matrix(cell)
        if m is None:
            return None
        try:
            value = m.contact(chrom, pos1, pos2, self.binsize)
        except LookupError:
            value = None
        cache[key] = value
        self._dirty.add(cell)
        return value

    def save(self) -> list[Path]:
        out = []
        for cell in sorted(self._dirty):
            p = self._path(cell)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(json.dumps(self._cache[cell], separators=(",", ":"), sort_keys=True))
            out.append(p)
        self._dirty.clear()
        return out

    def provenance(self) -> dict[str, Any]:
        reads = {
            cell: {"range_requests": m.reader.requests, "bytes_read": m.reader.bytes}
            for cell, m in self._open.items()
            if m is not None
        }
        return {
            "evidence": EVIDENCE,
            "matrices": {cell: download_url(a) for cell, a in self.matrices.items()},
            "binsize": self.binsize,
            "norm_preference": list(NORM_PREFERENCE),
            "reads": reads,  # empty when every answer came from the cache: no matrix was opened
            "cached_cells": {cell: len(self._load(cell)) for cell in self.matrices},
            "cache_bytes": {
                cell: (self._path(cell).stat().st_size if self._path(cell).exists() else 0)
                for cell in self.matrices
            },
        }
