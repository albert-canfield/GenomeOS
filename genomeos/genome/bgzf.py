"""Blocked gzip (BGZF): one compressed file that is still randomly accessible.

A BGZF file is an ordinary gzip file -- every reader that calls `gzip.open` on
it keeps working -- but it is written as a chain of independent members, each
holding at most 64 KiB of uncompressed data and carrying its own compressed
length in the gzip extra field. Because each member inflates on its own, a
companion index of (compressed offset, uncompressed offset) pairs turns any
uncompressed position into a seek plus one inflate. That index is the `.gzi`
file, written in the same layout htslib uses.

This is the format BAM has always used, and the project already decodes it from
a remote BAM in `bam_range.py`; here it is written and seeked locally, with the
standard library alone. The point is storage: a chromosome no longer has to be
kept twice, once compressed for streaming and once flat for random access.

Nothing here is lossy. The bytes that come back out of `BgzfReader` are the
bytes that went in, and `verify` says so by hash.
"""

from __future__ import annotations

import bisect
import hashlib
import os
import struct
import zlib
from collections.abc import Iterator
from pathlib import Path

# htslib writes at most 0xff00 uncompressed bytes per block, leaving room for
# the member to stay under the 65536-byte ceiling that BSIZE can express even
# when the data is incompressible.
BLOCK_SIZE = 0xFF00
MAX_MEMBER = 0x10000

# The empty member every BGZF file ends with, so a truncated file is detectable.
EOF_MARKER = (
    b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00"
    b"\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
)

_XFL = 0
_OS = 0xFF  # unknown, as htslib writes


def _member(payload: bytes, level: int) -> bytes:
    """One BGZF member: gzip header with the BC extra field, raw deflate, CRC and size."""
    co = zlib.compressobj(level, zlib.DEFLATED, -zlib.MAX_WBITS)
    body = co.compress(payload) + co.flush()
    total = len(body) + 26  # 12 header + 6 extra + 8 trailer
    if total > MAX_MEMBER:
        raise ValueError(f"BGZF member of {total} bytes exceeds {MAX_MEMBER}")
    head = struct.pack("<BBBBIBBHBBHH", 0x1F, 0x8B, 8, 4, 0, _XFL, _OS, 6, 0x42, 0x43, 2, total - 1)
    return head + body + struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload))


def _block_of(payload: bytes, level: int) -> tuple[bytes, int]:
    """A member for as much of `payload` as fits; returns it and how much was used."""
    take = len(payload)
    while True:
        try:
            return _member(payload[:take], level), take
        except ValueError:
            take = take * 3 // 4
            if take < 1024:  # pragma: no cover - incompressible beyond any sane input
                raise


def compress_file(src: str | Path, dest: str | Path, level: int = 6, progress=None) -> Path:
    """Write `src` as BGZF at `dest` and its `.gzi` index beside it.

    The destination is written to a `.part` file and renamed only once the whole
    input has been read, so an interrupted conversion never leaves a half file
    where a reader would find it.
    """
    src, dest = Path(src), Path(dest)
    tmp = dest.with_name(dest.name + ".part")
    gzi_tmp = dest.with_name(dest.name + ".gzi.part")
    entries: list[tuple[int, int]] = []  # (compressed offset, uncompressed offset) per block
    cpos = upos = 0
    with open(src, "rb") as fh, open(tmp, "wb") as out:
        while True:
            payload = fh.read(BLOCK_SIZE)
            if not payload:
                break
            member, used = _block_of(payload, level)
            if used < len(payload):  # pragma: no cover - only for incompressible input
                fh.seek(upos + used)
            out.write(member)
            entries.append((cpos, upos))
            cpos += len(member)
            upos += used
            if progress and len(entries) % 2048 == 0:
                progress(f"{upos / 1e6:.0f} MB in, {cpos / 1e6:.0f} MB out")
        out.write(EOF_MARKER)
    write_gzi(gzi_tmp, entries)
    gzi_tmp.rename(dest.with_name(dest.name + ".gzi"))
    tmp.rename(dest)
    return dest


def write_gzi(path: str | Path, entries: list[tuple[int, int]]) -> Path:
    """The htslib `.gzi` layout: a count, then one pair per block after the first."""
    path = Path(path)
    with open(path, "wb") as fh:
        fh.write(struct.pack("<Q", max(len(entries) - 1, 0)))
        for cpos, upos in entries[1:]:
            fh.write(struct.pack("<QQ", cpos, upos))
    return path


def read_gzi(path: str | Path) -> tuple[list[int], list[int]]:
    """Compressed and uncompressed block offsets, with the implicit first block restored."""
    raw = Path(path).read_bytes()
    (n,) = struct.unpack("<Q", raw[:8])
    if len(raw) < 8 + 16 * n:
        raise ValueError(f"{path}: gzi claims {n} entries but holds {(len(raw) - 8) // 16}")
    coffsets = [0]
    uoffsets = [0]
    for i in range(n):
        c, u = struct.unpack("<QQ", raw[8 + 16 * i : 24 + 16 * i])
        coffsets.append(c)
        uoffsets.append(u)
    return coffsets, uoffsets


def is_bgzf(path: str | Path) -> bool:
    """True when the file's first member carries the BC extra field BGZF is defined by."""
    path = Path(path)
    try:
        with open(path, "rb") as fh:
            head = fh.read(18)
    except OSError:
        return False
    if len(head) < 18 or head[:3] != b"\x1f\x8b\x08" or not head[3] & 4:
        return False
    xlen = struct.unpack("<H", head[10:12])[0]
    return xlen >= 6 and head[12:14] == b"BC" and struct.unpack("<H", head[14:16])[0] == 2


def gzi_path(path: str | Path) -> Path:
    p = Path(path)
    return p.with_name(p.name + ".gzi")


def index_file(path: str | Path) -> Path:
    """Build the `.gzi` for an existing BGZF file by walking its members."""
    path = Path(path)
    entries: list[tuple[int, int]] = []
    cpos = upos = 0
    size = path.stat().st_size
    with open(path, "rb") as fh:
        while cpos < size:
            head = fh.read(18)
            if len(head) < 18 or head[:2] != b"\x1f\x8b":
                break
            bsize = struct.unpack("<H", head[16:18])[0] + 1
            fh.seek(cpos + bsize - 4)
            (isize,) = struct.unpack("<I", fh.read(4))
            if isize:
                entries.append((cpos, upos))
            cpos += bsize
            upos += isize
            fh.seek(cpos)
    return write_gzi(gzi_path(path), entries)


class BgzfReader:
    """Random access into a BGZF file by uncompressed offset.

    `seek`/`read` behave like a plain binary file over the decompressed stream.
    One inflated block is kept, plus a small cache of recent ones, because the
    callers here read many short loci that fall in the same neighbourhood.
    """

    __slots__ = (
        "path",
        "_fh",
        "_coffsets",
        "_uoffsets",
        "_size",
        "_cache",
        "_order",
        "_cap",
        "_pos",
    )

    def __init__(self, path: str | Path, cache_blocks: int = 16) -> None:
        self.path = Path(path)
        gzi = gzi_path(self.path)
        if not gzi.exists():
            index_file(self.path)
        self._coffsets, self._uoffsets = read_gzi(gzi)
        self._fh = open(self.path, "rb")  # noqa: SIM115  (long-lived, closed by close())
        self._size = self._uncompressed_size()
        self._cache: dict[int, bytes] = {}
        self._order: list[int] = []
        self._cap = max(cache_blocks, 1)
        self._pos = 0

    def _uncompressed_size(self) -> int:
        """Total decompressed length: the last block's uncompressed offset plus its ISIZE."""
        last_c = self._coffsets[-1]
        self._fh.seek(last_c + 16)
        (bsize,) = struct.unpack("<H", self._fh.read(2))
        self._fh.seek(last_c + bsize + 1 - 4)
        (isize,) = struct.unpack("<I", self._fh.read(4))
        return self._uoffsets[-1] + isize

    def __len__(self) -> int:
        return self._size

    def close(self) -> None:
        self._fh.close()
        self._cache.clear()

    def __enter__(self) -> BgzfReader:
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _block(self, i: int) -> bytes:
        hit = self._cache.get(i)
        if hit is not None:
            return hit
        start = self._coffsets[i]
        self._fh.seek(start)
        head = self._fh.read(18)
        xlen = struct.unpack("<H", head[10:12])[0]
        bsize = struct.unpack("<H", head[16:18])[0] + 1
        self._fh.seek(start + 12 + xlen)
        body = self._fh.read(bsize - (12 + xlen) - 8)
        data = zlib.decompress(body, -zlib.MAX_WBITS)
        self._cache[i] = data
        self._order.append(i)
        if len(self._order) > self._cap:
            self._cache.pop(self._order.pop(0), None)
        return data

    def seek(self, offset: int, whence: int = os.SEEK_SET) -> int:
        if whence == os.SEEK_CUR:
            offset += self._pos
        elif whence == os.SEEK_END:
            offset += self._size
        self._pos = max(offset, 0)
        return self._pos

    def tell(self) -> int:
        return self._pos

    def read(self, n: int = -1) -> bytes:
        if n < 0:
            n = self._size - self._pos
        end = min(self._pos + n, self._size)
        if end <= self._pos:
            return b""
        out = self.read_range(self._pos, end)
        self._pos = end
        return out

    def read_range(self, start: int, end: int) -> bytes:
        """Bytes [start, end) of the decompressed stream, without moving the cursor."""
        if end <= start:
            return b""
        end = min(end, self._size)
        i = bisect.bisect_right(self._uoffsets, start) - 1
        pieces: list[bytes] = []
        pos = start
        while pos < end and i < len(self._uoffsets):
            block = self._block(i)
            base = self._uoffsets[i]
            a = pos - base
            b = min(end - base, len(block))
            if b > a:
                pieces.append(block[a:b])
                pos = base + b
            i += 1
        return b"".join(pieces)

    def blocks(self) -> Iterator[bytes]:
        """Every block in order, for a streaming pass that never holds the file."""
        for i in range(len(self._coffsets)):
            yield self._block(i)


def sha256_of(path: str | Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                return h.hexdigest()
            h.update(b)


def sha256_of_bgzf(path: str | Path, chunk: int = 1 << 24) -> str:
    """SHA-256 of what a BGZF file decompresses to, read through the seeking reader."""
    h = hashlib.sha256()
    with BgzfReader(path) as r:
        pos = 0
        total = len(r)
        while pos < total:
            h.update(r.read_range(pos, min(pos + chunk, total)))
            pos += chunk
    return h.hexdigest()


def verify(original: str | Path, bgzf: str | Path) -> tuple[bool, str, str]:
    """Compare the original file's SHA-256 with the BGZF file's decompressed SHA-256."""
    a = sha256_of(original)
    b = sha256_of_bgzf(bgzf)
    return a == b, a, b
