# SPDX-License-Identifier: AGPL-3.0-or-later
"""The shared range reader over kept-alive connections: one connection reused per host, a connection
the server dropped reopened without a failure, a failing host replaced by UCSC's mirror, sizes that
must agree, redirects followed and remembered, and a whole bigWig read identically either way."""

import http.client
import struct
import zlib

import pytest

from genomeos.attribution import bigwig as bw

FILE = bytes(range(256)) * 40  # the "remote file": 10,240 bytes


class FakeResponse:
    def __init__(self, status, body=b"", headers=None, will_close=False):
        self.status = status
        self._body = body
        self._headers = headers or {}
        self.will_close = will_close

    def read(self):
        return self._body

    def getheader(self, name, default=None):
        return self._headers.get(name, default)


class Server:
    """A scripted set of hosts: each host serves FILE, or fails, or redirects, as the plan says."""

    def __init__(self, plan=None, files=None):
        self.plan = plan or {}
        self.files = files or {}
        self.opened: list[str] = []
        self.requests: list[tuple[str, str]] = []

    def connection(self, host, timeout=None):
        self.opened.append(host)
        return FakeConnection(self, host)


class FakeConnection:
    def __init__(self, server, host):
        self.server, self.host = server, host
        self.pending = None
        self.closed = False
        self.used = 0

    def request(self, method, path, headers=None):
        self.server.requests.append((self.host, path))
        action = self.server.plan.get(self.host)
        if callable(action):
            action = action(self)
        self.used += 1
        if action == "fail":
            raise TimeoutError("timed out")
        if action == "drop" or self.closed:
            raise http.client.RemoteDisconnected("Remote end closed connection without response")
        if isinstance(action, tuple) and action[0] == "redirect":
            self.pending = FakeResponse(302, headers={"Location": action[1]})
            return
        if isinstance(action, int):
            self.pending = FakeResponse(action)
            return
        data = self.server.files.get(self.host, FILE)
        a, b = (int(x) for x in headers["Range"].split("=")[1].split("-"))
        body = data[a : b + 1]
        self.pending = FakeResponse(
            206, body, {"Content-Range": f"bytes {a}-{a + len(body) - 1}/{len(data)}"}
        )

    def getresponse(self):
        r, self.pending = self.pending, None
        return r

    def close(self):
        self.closed = True


@pytest.fixture
def server(monkeypatch):
    s = Server()
    monkeypatch.setattr(bw.http.client, "HTTPSConnection", s.connection)
    monkeypatch.setattr(bw.http.client, "HTTPConnection", s.connection)
    monkeypatch.setattr(bw.time, "sleep", lambda seconds: None)
    monkeypatch.setattr(bw, "_PROXIES", {})
    bw._POOL.conns = {}
    yield s
    bw._POOL.conns = {}


URL = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/track.bw"


def test_one_connection_is_reused_for_every_range(server):
    src = bw._Source(URL)
    for k in range(20):
        assert src.read(100 * k, 50) == FILE[100 * k : 100 * k + 50]
    assert server.opened == ["hgdownload.soe.ucsc.edu"]
    assert src.requests == 20 and src.bytes_fetched == 1000 and src.connections == 1
    # a second source on the same host, in the same thread, reuses the same connection
    other = bw._Source(URL)
    assert other.read(0, 10) == FILE[:10] and other.connections == 0 and len(server.opened) == 1


def test_a_dropped_keepalive_connection_is_reopened_without_a_retry(server):
    src = bw._Source(URL)
    assert src.read(0, 10) == FILE[:10]
    calls = {"n": 0}

    def drop_once(conn):
        calls["n"] += 1
        return "drop" if calls["n"] == 1 else None

    server.plan["hgdownload.soe.ucsc.edu"] = drop_once
    assert src.read(5000, 64) == FILE[5000:5064]
    assert server.opened == ["hgdownload.soe.ucsc.edu", "hgdownload.soe.ucsc.edu"]
    assert src.mirror_requests == 0  # reopened on the same host, not counted as a failure


def test_a_connection_dropped_mid_stream_is_never_reused(server, monkeypatch):
    src = bw._Source(URL)
    assert src.read(0, 10) == FILE[:10]
    conn = bw._POOL.conns[("https", "hgdownload.soe.ucsc.edu")]

    def broken_read():
        raise http.client.IncompleteRead(b"partial")

    real = conn.getresponse

    def half(*args):
        r = real()
        r.read = broken_read
        return r

    monkeypatch.setattr(conn, "getresponse", half)
    # the half-read connection is discarded; the retry goes to the mirror on a new connection
    assert src.read(20, 10) == FILE[20:30]
    assert ("https", "hgdownload.soe.ucsc.edu") not in bw._POOL.conns or bw._POOL.conns[
        ("https", "hgdownload.soe.ucsc.edu")
    ] is not conn
    assert src.mirror_requests == 1


def test_a_failing_host_fails_over_to_the_mirror(server):
    server.plan["hgdownload.soe.ucsc.edu"] = "fail"
    src = bw._Source(URL)
    assert src.read(64, 32) == FILE[64:96]
    assert src.mirror_requests == 1 and src.total == len(FILE)
    assert server.requests[-1] == ("hgdownload2.soe.ucsc.edu", "/gbdb/hg38/track.bw")


def test_hosts_that_disagree_on_the_file_raise(server):
    src = bw._Source(URL)
    src.read(0, 8)
    server.plan["hgdownload.soe.ucsc.edu"] = "fail"
    server.files["hgdownload2.soe.ucsc.edu"] = FILE + b"extra"
    with pytest.raises(ValueError, match="disagree"):
        src.read(0, 8)


def test_server_errors_are_retried_and_then_raised(server):
    server.plan["hgdownload.soe.ucsc.edu"] = 503
    server.plan["hgdownload2.soe.ucsc.edu"] = 503
    src = bw._Source(URL, retries=3)
    with pytest.raises(OSError, match="server error"):
        src.read(0, 8)
    assert len(server.requests) == 3


def test_redirects_are_followed_and_remembered(server):
    server.plan["www.encodeproject.org"] = ("redirect", "https://bucket.example.org/file.bigWig?sig=1")
    src = bw._Source("https://www.encodeproject.org/files/X/@@download/X.bigWig")
    assert src.read(10, 5) == FILE[10:15]
    assert src.read(20, 5) == FILE[20:25]
    hosts = [h for h, _ in server.requests]
    assert hosts == ["www.encodeproject.org", "bucket.example.org", "bucket.example.org"]
    # an expired target is resolved again from the original address
    calls = {"n": 0}

    def expire_once(conn):
        calls["n"] += 1
        return 403 if calls["n"] == 1 else None

    server.plan["bucket.example.org"] = expire_once
    assert src.read(30, 5) == FILE[30:35]
    assert [h for h, _ in server.requests][-3:] == [
        "bucket.example.org",
        "www.encodeproject.org",
        "bucket.example.org",
    ]


def _bigwig_file() -> bytes:
    """A one-chromosome bigWig small enough to build by hand: a fixedStep section of 100 values."""
    values = [float(i % 7) for i in range(100)]
    section = struct.pack("<IIIIIBBH", 0, 1000, 1100, 1, 1, bw.FIXED_STEP, 0, 100) + struct.pack(
        "<100f", *values
    )
    data = zlib.compress(section)
    header_size, chrom_tree_at = 64, 64
    key = b"chr1"
    chrom_tree = (
        struct.pack("<IIIIQQ", bw.CHROM_TREE_MAGIC, 1, 4, 8, 1, 0)
        + struct.pack("<BBH", 1, 0, 1)
        + key
        + struct.pack("<II", 0, 5000)
    )
    data_at = chrom_tree_at + len(chrom_tree)
    index_at = data_at + len(data)
    rtree = struct.pack("<IIQIIIIQII", bw.RTREE_MAGIC, 1, 1, 0, 1000, 0, 1100, index_at, 1, 0)
    leaf = struct.pack("<BBH", 1, 0, 1) + struct.pack("<IIIIQQ", 0, 1000, 0, 1100, data_at, len(data))
    leaf += b"\0" * 64
    header = struct.pack(
        "<IHHQQQHHQQIQ", bw.BIGWIG_MAGIC, 4, 0, chrom_tree_at, data_at, index_at, 0, 0, 0, 0, 1 << 15, 0
    )
    assert len(header) == header_size
    return header + chrom_tree + data + rtree + leaf


def test_a_whole_bigwig_reads_the_same_from_disk_and_over_kept_alive_ranges(server, tmp_path):
    blob = _bigwig_file()
    path = tmp_path / "t.bw"
    path.write_bytes(blob)
    server.files["hgdownload.soe.ucsc.edu"] = blob
    local = bw.BigWig(path).summarise("chr1", [(1000, 1050), (1090, 1100)], 3.0)
    remote_reader = bw.BigWig(URL)
    remote = remote_reader.summarise("chr1", [(1000, 1050), (1090, 1100)], 3.0)
    assert [s.as_dict() for s in local] == [s.as_dict() for s in remote]
    assert local[0].bases == 50 and remote_reader.src.connections == 1
