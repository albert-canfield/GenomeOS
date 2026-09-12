# SPDX-License-Identifier: AGPL-3.0-or-later
"""A coordinate-sorted, indexed BAM read by HTTP ranges, standard library only.

A 300x whole-genome BAM is 600 GB. Its index (.bai, 12 MB) says where each 16 kb window of each
chromosome starts in the compressed file and, in its two pseudo-bins, how many reads mapped to each
reference and how many never mapped. With that, the reads over one window, or a sample of the unmapped
tail, cost a few megabytes of range requests and nothing is downloaded whole. BGZF blocks are gzip
members with the block size in an extra field, so the decoder is the standard library's.

Read once, kept nowhere: the caller keeps counts. Used by the telomere estimate; general enough for
any window of any indexed BAM that answers range requests.
"""

from __future__ import annotations

import gzip
import struct
import urllib.request
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

USER_AGENT = "GenomeOS/0.1 (stream)"
SEQ = "=ACMGRSVTWYHKDBN"
PSEUDO_BIN = 37450
INTERVAL = 16_384
EOF_BLOCK = 28  # the empty BGZF block that ends every BAM


@dataclass
class RefIndex:
    name: str
    length: int
    mapped: int | None = None
    unmapped: int | None = None
    intervals: tuple[int, ...] = ()
    chunk_end: int = 0  # the largest compressed offset any of this reference's chunks reach


@dataclass
class BamIndex:
    refs: list[RefIndex]
    no_coor: int | None = None  # reads with no coordinate at all (the unmapped tail)
    by_name: dict[str, int] = field(default_factory=dict)

    @property
    def mapped_total(self) -> int | None:
        counts = [r.mapped for r in self.refs if r.mapped is not None]
        return sum(counts) if counts else None

    @property
    def tail_start(self) -> int:
        """Compressed offset where the unmapped tail begins (after the last mapped chunk)."""
        return max(r.chunk_end for r in self.refs) >> 16


def parse_header(blob: bytes) -> tuple[str, list[tuple[str, int]]]:
    """The SAM text and the references (name, length) from the decompressed head of a BAM."""
    if blob[:4] != b"BAM\x01":
        raise ValueError("not a BAM")
    (l_text,) = struct.unpack("<i", blob[4:8])
    text = blob[8 : 8 + l_text].decode(errors="replace")
    p = 8 + l_text
    (n_ref,) = struct.unpack("<i", blob[p : p + 4])
    p += 4
    refs = []
    for _ in range(n_ref):
        (ln,) = struct.unpack("<i", blob[p : p + 4])
        p += 4
        name = blob[p : p + ln - 1].decode()
        p += ln
        (length,) = struct.unpack("<i", blob[p : p + 4])
        p += 4
        refs.append((name, length))
    return text, refs


def parse_bai(blob: bytes, refs: list[tuple[str, int]] | None = None) -> BamIndex:
    """The linear index, the per-reference read counts (pseudo-bin 37450) and the unmapped count."""
    if blob[:4] != b"BAI\x01":
        raise ValueError("not a BAI")
    (n_ref,) = struct.unpack("<i", blob[4:8])
    p = 8
    out: list[RefIndex] = []
    for r in range(n_ref):
        name, length = refs[r] if refs and r < len(refs) else (str(r), 0)
        ri = RefIndex(name, length)
        (n_bin,) = struct.unpack("<i", blob[p : p + 4])
        p += 4
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack("<Ii", blob[p : p + 8])
            p += 8
            chunks = [struct.unpack("<QQ", blob[p + 16 * i : p + 16 * i + 16]) for i in range(n_chunk)]
            p += 16 * n_chunk
            if bin_id == PSEUDO_BIN and len(chunks) == 2:
                ri.mapped, ri.unmapped = chunks[1]
            else:
                for _, end in chunks:
                    ri.chunk_end = max(ri.chunk_end, end)
        (n_intv,) = struct.unpack("<i", blob[p : p + 4])
        p += 4
        ri.intervals = struct.unpack(f"<{n_intv}Q", blob[p : p + 8 * n_intv])
        p += 8 * n_intv
        out.append(ri)
    no_coor = struct.unpack("<Q", blob[p : p + 8])[0] if p + 8 <= len(blob) else None
    return BamIndex(out, no_coor, {r.name: i for i, r in enumerate(out)})


def decode_bgzf(raw: bytes) -> tuple[bytes, int]:
    """Every complete BGZF block in `raw`, decompressed; also how many compressed bytes were used."""
    out = bytearray()
    off = 0
    while off + 18 <= len(raw) and raw[off : off + 2] == b"\x1f\x8b":
        bsize = struct.unpack("<H", raw[off + 16 : off + 18])[0] + 1
        if off + bsize > len(raw):
            break
        out += gzip.decompress(raw[off : off + bsize])
        off += bsize
    return bytes(out), off


@dataclass(slots=True)
class Read:
    ref_id: int
    pos: int
    flag: int
    seq: str


def records(data: bytes) -> Iterator[Read]:
    """BAM alignment records from decompressed bytes that start at a record boundary."""
    q = 0
    n = len(data)
    while q + 36 <= n:
        (bs,) = struct.unpack("<i", data[q : q + 4])
        rec = data[q + 4 : q + 4 + bs]
        if len(rec) < bs or bs < 32:
            return
        q += 4 + bs
        ref_id, pos = struct.unpack("<ii", rec[0:8])
        l_rn = rec[8]
        n_cig = struct.unpack("<H", rec[12:14])[0]
        flag = struct.unpack("<H", rec[14:16])[0]
        (l_seq,) = struct.unpack("<i", rec[16:20])
        o = 32 + l_rn + 4 * n_cig
        packed = rec[o : o + (l_seq + 1) // 2]
        seq = "".join(SEQ[x >> 4] + SEQ[x & 15] for x in packed)[:l_seq]
        yield Read(ref_id, pos, flag, seq)


class RemoteBam:
    """One indexed BAM behind a URL (or a local path), read by ranges."""

    def __init__(self, source: str, index: bytes | str | Path, head_bytes: int = 1 << 20) -> None:
        self.source = source
        self.remote = source.startswith(("http://", "https://"))
        self.bytes_fetched = 0
        self.requests = 0
        head, _ = decode_bgzf(self.read(0, head_bytes))
        self.header_text, refs = parse_header(head)
        blob = index if isinstance(index, bytes) else Path(index).read_bytes()
        self.index = parse_bai(blob, refs)
        self.size = self._size()

    def _size(self) -> int:
        if not self.remote:
            return Path(self.source).stat().st_size
        req = urllib.request.Request(self.source, method="HEAD", headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310
            return int(r.headers.get("Content-Length", "0"))

    def read(self, offset: int, size: int) -> bytes:
        self.requests += 1
        self.bytes_fetched += size
        if not self.remote:
            with open(self.source, "rb") as fh:
                fh.seek(offset)
                return fh.read(size)
        req = urllib.request.Request(
            self.source,
            headers={"Range": f"bytes={offset}-{offset + size - 1}", "User-Agent": USER_AGENT},
        )
        with urllib.request.urlopen(req, timeout=900) as r:  # noqa: S310
            return r.read()

    def reads(self, chrom: str, start: int, end: int) -> Iterator[Read]:
        """Reads whose leftmost position falls in [start, end) on `chrom` (0-based)."""
        ref = self.index.by_name[chrom]
        ri = self.index.refs[ref]
        first = next((v for v in ri.intervals[start // INTERVAL :] if v), None)
        if first is None:
            return
        later = [v for v in ri.intervals[end // INTERVAL + 1 :] if v]
        for nxt in self.index.refs[ref + 1 :]:
            later += [v for v in nxt.intervals if v][:1]
            if later:
                break
        c0 = first >> 16
        c1 = (min(later) >> 16) if later else self.index.tail_start
        raw = self.read(c0, max(c1 - c0, 1) + 65_536)
        data, _ = decode_bgzf(raw)
        for rd in records(data[first & 0xFFFF :]):
            if rd.ref_id != ref:
                if rd.ref_id > ref or rd.ref_id < 0:
                    return
                continue
            if rd.pos >= end:
                return
            if rd.pos >= start:
                yield rd

    def unmapped_sample(self, max_bytes: int) -> Iterator[Read]:
        """The first `max_bytes` of the unmapped tail (reads with no coordinate), decoded."""
        start = max(r.chunk_end for r in self.index.refs)
        c0 = start >> 16
        raw = self.read(c0, min(max_bytes, max(1, self.size - EOF_BLOCK - c0)))
        data, _ = decode_bgzf(raw)
        yield from records(data[start & 0xFFFF :])
