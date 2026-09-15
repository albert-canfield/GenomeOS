# SPDX-License-Identifier: AGPL-3.0-or-later
"""A bigWig reader over HTTP range requests, standard library only.

The Zoonomia phyloP track for hg38 is a 9.6 GB bigWig. Nothing here downloads it:
the chromosome tree and the R-tree index are read with range requests, the data
sections that overlap the intervals asked for are fetched in a few coalesced
ranges, decompressed and summarised, and the bytes are discarded. Asking for the
UNKNOWN blocks of chromosome 21 moves about 40 MB.

Format reference: Kent et al. 2010, "BigWig and BigBed", and the UCSC kent source
(bbiFile.h, bwgInternal.h). Little-endian files only, which is what UCSC writes.
"""

from __future__ import annotations

import http.client
import struct
import sys
import threading
import time
import urllib.parse
import urllib.request
import zlib
from array import array
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path
from typing import Any

BIGWIG_MAGIC = 0x888FFC26
CHROM_TREE_MAGIC = 0x78CA8C91
RTREE_MAGIC = 0x2468ACE0
USER_AGENT = "GenomeOS/0.9 (bigwig range reader; one chromosome's blocks at a time)"

# section types inside a bigWig data block
BEDGRAPH, VARIABLE_STEP, FIXED_STEP = 1, 2, 3


# UCSC serves the same files from a second download host; a request that fails on one is retried on
# the other (genomeos-h1 checked the sizes on 2026-09-14, and every response's total size must agree).
MIRRORS = {
    "hgdownload.soe.ucsc.edu": "hgdownload2.soe.ucsc.edu",
    "hgdownload2.soe.ucsc.edu": "hgdownload.soe.ucsc.edu",
}
MAX_REDIRECTS = 5
MAX_DELAY = 30.0
_POOL = threading.local()  # (scheme, host) -> an open connection, one set per thread
_DROPPED = (
    http.client.RemoteDisconnected,
    ConnectionResetError,
    BrokenPipeError,
    http.client.CannotSendRequest,
)


_PROXIES: dict[str, str] | None = None


def _proxied(scheme: str) -> bool:
    """A proxy in the environment: requests are left to urllib, which knows how to use one."""
    global _PROXIES
    if _PROXIES is None:
        _PROXIES = urllib.request.getproxies()
    return bool(_PROXIES.get(scheme))


class _Source:
    """Random access to a local file or a URL, by (offset, size).

    Over HTTP every range used to open a new TLS connection. genomeos-h1 profiled their panel run
    (2026-09-14): the handshake costs 0.5 s when UCSC is quiet and 4 to 17 s when it is loaded, and
    one connection kept open per host took a run from 31,099 s to 1,392 s with identical numbers.
    The same approach is used here (their `KeepAliveSource` in attribution/human_panel.py): one
    persistent connection per host and thread, reused for every range; a dropped keep-alive
    connection is reopened at once without counting as a failure; any other failure is retried
    with backoff, alternating with UCSC's second download host, and the file's total size must be
    the same on every response. Redirects (ENCODE's download links) are followed and the final
    address remembered. Interface and results are unchanged: `read`, `close`, `requests`,
    `bytes_fetched`, and the same bytes for the same (offset, size).
    """

    def __init__(self, source: str | Path, timeout: int = 180, retries: int = 6):
        self.remote = isinstance(source, str) and source.startswith(("http://", "https://"))
        self.source = str(source)
        self.timeout = timeout
        self.retries = retries
        self.bytes_fetched = 0
        self.requests = 0
        self.connections = 0  # connections this source opened (a pooled one reused costs none)
        self.mirror_requests = 0
        self.total: int | None = None
        self._resolved: str | None = None  # where a redirect led, reused until it fails
        self._fh = None if self.remote else open(self.source, "rb")  # noqa: SIM115

    # -- connections ----------------------------------------------------------------------------
    def _conn(self, scheme: str, host: str, fresh: bool = False) -> tuple[http.client.HTTPConnection, bool]:
        """The pooled connection to a host, and whether it was reused rather than opened now."""
        pool = getattr(_POOL, "conns", None)
        if pool is None:
            pool = _POOL.conns = {}
        key = (scheme, host)
        conn = pool.get(key)
        if conn is not None and not fresh:
            return conn, True
        if conn is not None:
            conn.close()
        cls = http.client.HTTPSConnection if scheme == "https" else http.client.HTTPConnection
        conn = pool[key] = cls(host, timeout=self.timeout)
        self.connections += 1
        return conn, False

    def _drop(self, scheme: str, host: str) -> None:
        pool = getattr(_POOL, "conns", None) or {}
        conn = pool.pop((scheme, host), None)
        if conn is not None:
            conn.close()

    def _once(self, url: str, headers: dict[str, str]) -> tuple[int, bytes, Any, str]:
        """One request, following redirects. A kept-alive connection the server closed while idle
        is reopened once on the spot; a connection that fails in any other way is discarded, never
        reused half-read. Returns the status, the body, the response and the final URL."""
        for _ in range(MAX_REDIRECTS + 1):
            u = urllib.parse.urlsplit(url)
            path = (u.path or "/") + (f"?{u.query}" if u.query else "")
            for again in (False, True):
                conn, reused = self._conn(u.scheme, u.netloc, fresh=again)
                try:
                    conn.request("GET", path, headers={"User-Agent": USER_AGENT, **headers})
                    r = conn.getresponse()
                    data = r.read()
                    break
                except _DROPPED:
                    self._drop(u.scheme, u.netloc)
                    if again or not reused:
                        raise
                except BaseException:
                    self._drop(u.scheme, u.netloc)
                    raise
            if r.will_close:
                self._drop(u.scheme, u.netloc)
            location = r.getheader("Location")
            if r.status in (301, 302, 303, 307, 308) and location:
                url = urllib.parse.urljoin(url, location)
                continue
            return r.status, data, r, url
        raise OSError(f"{self.source}: more than {MAX_REDIRECTS} redirects")

    def _read_urllib(self, offset: int, size: int) -> bytes:  # pragma: no cover - proxied networks
        req = urllib.request.Request(
            self.source,
            headers={"Range": f"bytes={offset}-{offset + size - 1}", "User-Agent": USER_AGENT},
        )
        with urllib.request.urlopen(req, timeout=self.timeout) as r:  # noqa: S310
            data = r.read()
        if len(data) < size and r.status != 206:
            raise OSError(f"server ignored the range request ({r.status})")
        return data

    def _read_remote(self, offset: int, size: int) -> bytes:
        headers = {"Range": f"bytes={offset}-{offset + size - 1}"}
        scheme = urllib.parse.urlsplit(self.source).scheme
        delay = 2.0
        for attempt in range(self.retries):
            try:
                if _proxied(scheme):  # pragma: no cover - proxied networks
                    return self._read_urllib(offset, size)
                url = self._resolved or self.source
                host = urllib.parse.urlsplit(url).netloc
                mirrored = attempt % 2 == 1 and host in MIRRORS
                if mirrored:
                    url = url.replace(f"//{host}/", f"//{MIRRORS[host]}/", 1)
                status, data, r, final = self._once(url, headers)
                if status >= 500:
                    raise OSError(f"server error {status}")
                if status in (401, 403, 404, 410) and self._resolved:
                    self._resolved = None  # a remembered redirect target expired: resolve again
                    raise OSError(f"redirect target refused ({status})")
                if status not in (200, 206):
                    raise OSError(f"range request failed ({status})")
                if status == 200 and (offset != 0 or len(data) != size):
                    raise OSError(f"server ignored the range request ({status})")
                total = (r.getheader("Content-Range") or "").rpartition("/")[2]
                if total.isdigit():
                    if self.total is None:
                        self.total = int(total)
                    elif int(total) != self.total:
                        raise ValueError(
                            f"{self.source}: hosts disagree on the file size ({total} against {self.total})"
                        )
                if final != url and not mirrored:
                    self._resolved = final
                self.mirror_requests += mirrored
                return data
            except (OSError, http.client.HTTPException):
                if attempt == self.retries - 1:
                    raise
                time.sleep(delay)
                delay = min(delay * 2, MAX_DELAY)
        raise OSError("unreachable")  # pragma: no cover

    def read(self, offset: int, size: int) -> bytes:
        self.requests += 1
        self.bytes_fetched += size
        if not self.remote:
            self._fh.seek(offset)
            return self._fh.read(size)
        return self._read_remote(offset, size)

    def close(self) -> None:
        if self._fh:
            self._fh.close()


@dataclass(frozen=True)
class Header:
    version: int
    zoom_levels: int
    chrom_tree_offset: int
    full_data_offset: int
    full_index_offset: int
    uncompress_buf_size: int


@dataclass(frozen=True)
class LeafItem:
    start_chrom: int
    start_base: int
    end_chrom: int
    end_base: int
    offset: int
    size: int


@dataclass
class IntervalStats:
    """What the track says over one interval: bases with a value, their sum, max, and how
    many reach the threshold."""

    start: int
    end: int
    bases: int = 0
    total: float = 0.0
    maximum: float = float("-inf")
    above: int = 0

    @property
    def mean(self) -> float | None:
        return self.total / self.bases if self.bases else None

    @property
    def fraction_above(self) -> float | None:
        return self.above / self.bases if self.bases else None

    def as_dict(self) -> dict:
        return {
            "bases": self.bases,
            "mean": round(self.mean, 4) if self.bases else None,
            "max": round(self.maximum, 3) if self.bases else None,
            "above": self.above,
            "fraction_above": round(self.fraction_above, 4) if self.bases else None,
        }


def decode_section(raw: bytes) -> list[tuple[int, int, array]]:
    """One decompressed data section as runs of (start, step, values).

    For fixedStep with span == step the run is dense: value i covers
    ``start + i*step .. + step``. bedGraph and variableStep sections are expanded to
    the same shape per item so callers have one code path; their runs carry a
    single value each (start, span, [v]).
    """
    chrom_start, chrom_end, step, span, kind, count = struct.unpack("<xxxxIIIIBxH", raw[:24])
    body = raw[24:]
    if kind == FIXED_STEP:
        vals = array("f")
        vals.frombytes(body[: 4 * count])
        if sys.byteorder == "big":  # pragma: no cover
            vals.byteswap()
        if span == step:
            return [(chrom_start, step, vals)]
        return [(chrom_start + i * step, span, array("f", [v])) for i, v in enumerate(vals)]
    out = []
    if kind == VARIABLE_STEP:
        for i in range(count):
            s, v = struct.unpack_from("<If", body, 8 * i)
            out.append((s, span, array("f", [v])))
        return out
    if kind == BEDGRAPH:
        for i in range(count):
            s, e, v = struct.unpack_from("<IIf", body, 12 * i)
            out.append((s, e - s, array("f", [v])))
        return out
    raise ValueError(f"unknown bigWig section type {kind}")


def coalesce(items: list[LeafItem], max_gap: int = 1 << 16, max_size: int = 8 << 20) -> list[list[LeafItem]]:
    """Group data blocks so one range request serves many: sorted by offset, joined while the
    gap between them stays under `max_gap` and the request under `max_size`."""
    groups: list[list[LeafItem]] = []
    for it in sorted(items, key=lambda i: i.offset):
        if groups:
            g = groups[-1]
            first, last = g[0], g[-1]
            if (
                it.offset - (last.offset + last.size) <= max_gap
                and it.offset + it.size - first.offset <= max_size
            ):
                g.append(it)
                continue
        groups.append([it])
    return groups


class BigWig:
    """Open a bigWig by path or URL; summarise intervals of one chromosome."""

    def __init__(self, source: str | Path):
        self.src = _Source(source)
        h = self.src.read(0, 64)
        magic, version, zoom, ct, fd, fi, _fc, _dfc, _asql, _ts, ubs, _res = struct.unpack("<IHHQQQHHQQIQ", h)
        if magic != BIGWIG_MAGIC:
            raise ValueError(f"not a little-endian bigWig: magic {magic:#x}")
        self.header = Header(version, zoom, ct, fd, fi, ubs)
        self.chroms: dict[str, tuple[int, int]] = {}
        self._read_chrom_tree()
        (_m, self._rtree_block_size, _n, _sc, _sb, _ec, _eb, _eo, self._items_per_slot, _r) = struct.unpack(
            "<IIQIIIIQII", self.src.read(fi, 48)
        )
        self._node_cache: dict[int, bytes] = {}

    def close(self) -> None:
        self.src.close()

    # -- the chromosome B+ tree -------------------------------------------------------------
    def _read_chrom_tree(self) -> None:
        off = self.header.chrom_tree_offset
        magic, block_size, key_size, _val_size, _count, _res = struct.unpack(
            "<IIIIQQ", self.src.read(off, 32)
        )
        if magic != CHROM_TREE_MAGIC:
            raise ValueError("bad chromosome tree")
        self._walk_chrom_node(off + 32, block_size, key_size)

    def _walk_chrom_node(self, off: int, block_size: int, key_size: int) -> None:
        item = key_size + 8
        raw = self.src.read(off, 4 + block_size * item)
        is_leaf, _r, count = struct.unpack("<BBH", raw[:4])
        for i in range(count):
            p = 4 + i * item
            key = raw[p : p + key_size].split(b"\0", 1)[0].decode("ascii")
            if is_leaf:
                cid, size = struct.unpack_from("<II", raw, p + key_size)
                self.chroms[key] = (cid, size)
            else:
                (child,) = struct.unpack_from("<Q", raw, p + key_size)
                self._walk_chrom_node(child, block_size, key_size)

    # -- the R-tree over data blocks ------------------------------------------------------------
    def _node(self, off: int) -> bytes:
        raw = self._node_cache.get(off)
        if raw is None:
            raw = self.src.read(off, 4 + self._rtree_block_size * 32)
            self._node_cache[off] = raw
        return raw

    def leaf_items(self, chrom: str, intervals: list[tuple[int, int]]) -> list[LeafItem]:
        """Every data block overlapping any of the (sorted, non-overlapping) intervals."""
        cid = self.chroms[chrom][0]
        starts = [s for s, _ in intervals]
        ends = [e for _, e in intervals]

        def overlaps(sc: int, sb: int, ec: int, eb: int) -> bool:
            if ec < cid or sc > cid:
                return False
            lo = sb if sc == cid else 0
            hi = eb if ec == cid else 1 << 62
            i = bisect_right(ends, lo)
            return i < len(starts) and starts[i] < hi

        out: list[LeafItem] = []

        def walk(off: int) -> None:
            raw = self._node(off)
            is_leaf, _r, count = struct.unpack("<BBH", raw[:4])
            if is_leaf:
                for i in range(count):
                    sc, sb, ec, eb, doff, dsize = struct.unpack_from("<IIIIQQ", raw, 4 + 32 * i)
                    if overlaps(sc, sb, ec, eb):
                        out.append(LeafItem(sc, sb, ec, eb, doff, dsize))
            else:
                for i in range(count):
                    sc, sb, ec, eb, child = struct.unpack_from("<IIIIQ", raw, 4 + 24 * i)
                    if overlaps(sc, sb, ec, eb):
                        walk(child)

        walk(self.header.full_index_offset + 48)
        return out

    # -- summaries ------------------------------------------------------------------------------
    def summarise(
        self,
        chrom: str,
        intervals: list[tuple[int, int]],
        threshold: float,
        progress=None,
    ) -> list[IntervalStats]:
        """Per interval: bases with a value, mean, max and the count at or above `threshold`.

        Intervals may be given in any order; they must not overlap each other.
        """
        order = sorted(range(len(intervals)), key=lambda i: intervals[i])
        ivs = [intervals[i] for i in order]
        stats = [IntervalStats(s, e) for s, e in ivs]
        if not ivs or chrom not in self.chroms:
            return [stats[order.index(i)] for i in range(len(intervals))] if ivs else []
        cid = self.chroms[chrom][0]
        ends = [e for _, e in ivs]
        items = self.leaf_items(chrom, ivs)
        groups = coalesce(items)
        done = 0
        for g in groups:
            base = g[0].offset
            blob = self.src.read(base, g[-1].offset + g[-1].size - base)
            for it in g:
                raw = blob[it.offset - base : it.offset - base + it.size]
                if self.header.uncompress_buf_size:
                    raw = zlib.decompress(raw)
                sec_cid = struct.unpack_from("<I", raw, 0)[0]
                if sec_cid != cid:
                    continue
                for start, step, vals in decode_section(raw):
                    self._accumulate(start, step, vals, ivs, ends, stats, threshold)
            done += len(g)
            if progress:
                progress(done, len(items), self.src.bytes_fetched)
        # back to the caller's order
        inv = [0] * len(order)
        for rank, i in enumerate(order):
            inv[i] = rank
        return [stats[inv[i]] for i in range(len(intervals))]

    @staticmethod
    def _accumulate(
        start: int,
        step: int,
        vals: array,
        ivs: list[tuple[int, int]],
        ends: list[int],
        stats: list[IntervalStats],
        threshold: float,
    ) -> None:
        n = len(vals)
        run_end = start + n * step
        i = bisect_right(ends, start)
        # the track stores float32; round the threshold the same way so 2.27 counts as 2.27
        ge = array("f", [threshold])[0].__le__
        while i < len(ivs) and ivs[i][0] < run_end:
            s, e = ivs[i]
            lo, hi = max(s, start), min(e, run_end)
            if lo < hi:
                if step == 1:
                    part = vals[lo - start : hi - start]
                else:
                    part = vals[(lo - start) // step : (hi - start + step - 1) // step]
                st = stats[i]
                covered = (hi - lo) if step == 1 else len(part) * step
                st.bases += covered
                st.total += sum(part) * (1 if step == 1 else step)
                m = max(part)
                if m > st.maximum:
                    st.maximum = m
                st.above += sum(map(ge, part)) * (1 if step == 1 else step)
            i += 1
