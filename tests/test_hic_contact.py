# SPDX-License-Identifier: AGPL-3.0-or-later
"""The .hic range reader on a synthetic matrix: header, index, block, balancing and expected value.

No network: a small version 8 .hic file is built in memory, handed to the reader through a fake range
reader, and the contacts it returns are checked against the numbers that were written in. The test
also checks that the reader touches only a fraction of the file, which is the whole point of it.
"""

import json
import math
import struct
import zlib

import pytest

from genomeos.genome import hic_contact as hc


def s(text: str) -> bytes:
    return text.encode() + b"\x00"


def block_bytes(records, bin_x_offset, bin_y_offset):
    """One v8 block, layout 1 (row lists) with short counts, as juicer writes it."""
    rows = {}
    for (x, y), v in sorted(records.items()):
        rows.setdefault(y, []).append((x, v))
    out = struct.pack("<i", len(records))
    out += struct.pack("<ii", bin_x_offset, bin_y_offset)
    out += bytes([0])  # short counts
    out += bytes([1])  # layout 1
    out += struct.pack("<h", len(rows))
    for y, cells in sorted(rows.items()):
        out += struct.pack("<h", y - bin_y_offset)
        out += struct.pack("<h", len(cells))
        for x, v in cells:
            out += struct.pack("<h", x - bin_x_offset) + struct.pack("<h", int(v))
    return zlib.compress(out)


BINSIZE = 5_000
BLOCK_BIN_COUNT = 100
BLOCK_COLUMN_COUNT = 10
CHROMS = [("ALL", 1000), ("1", 2_000_000)]
RECORDS = {(0, 0): 500, (0, 4): 40, (4, 4): 300, (10, 60): 7}
WEIGHTS = [1.0] * 400
EXPECTED = [100.0, 50.0, 25.0, 12.0, 6.0] + [1.0] * 395
SCALE_FACTOR = 2.0


def synthetic_hic() -> bytes:
    """A version 8 .hic file with one chromosome, one resolution, one block and a KR vector."""
    WEIGHTS[0], WEIGHTS[4], WEIGHTS[10] = 2.0, 0.5, 4.0
    WEIGHTS[60] = float("nan")  # a bin the balancing left undefined
    body = block_bytes(RECORDS, 0, 0)
    parts: list[bytes] = []
    header = b"HIC\x00" + struct.pack("<i", 8)
    header += struct.pack("<q", 0)  # footer position, filled in below
    header += s("synthetic")
    header += struct.pack("<i", 0)  # no attributes
    header += struct.pack("<i", len(CHROMS))
    for name, size in CHROMS:
        header += s(name) + struct.pack("<i", size)
    header += struct.pack("<i", 1) + struct.pack("<i", BINSIZE)  # one bp resolution
    header += struct.pack("<i", 0)  # no fragment resolutions
    block_position = len(header)
    parts.append(header)
    parts.append(body)

    matrix_position = block_position + len(body)
    matrix = struct.pack("<iii", 1, 1, 1)
    matrix += s("BP") + struct.pack("<i", 0)
    matrix += struct.pack("<ffff", 0.0, 0.0, 0.0, 0.0)
    matrix += struct.pack("<iiii", BINSIZE, BLOCK_BIN_COUNT, BLOCK_COLUMN_COUNT, 1)
    matrix += struct.pack("<i", 0) + struct.pack("<q", block_position) + struct.pack("<i", len(body))
    parts.append(matrix)

    norm_position = matrix_position + len(matrix)
    norm = struct.pack("<i", len(WEIGHTS)) + struct.pack(f"<{len(WEIGHTS)}d", *WEIGHTS)
    parts.append(norm)

    footer_position = norm_position + len(norm)
    footer = struct.pack("<i", 0)  # bytes of the v5 footer, unused
    footer += struct.pack("<i", 1)
    footer += s("1_1") + struct.pack("<q", matrix_position) + struct.pack("<i", len(matrix))
    footer += struct.pack("<i", 0)  # no unnormalised expected value vectors
    footer += struct.pack("<i", 1)  # one normalised expected value vector
    footer += s("KR") + s("BP") + struct.pack("<i", BINSIZE)
    footer += struct.pack("<i", len(EXPECTED)) + struct.pack(f"<{len(EXPECTED)}d", *EXPECTED)
    footer += struct.pack("<i", 1) + struct.pack("<i", 1) + struct.pack("<d", SCALE_FACTOR)
    footer += struct.pack("<i", 1)  # one normalisation vector
    footer += s("KR") + struct.pack("<i", 1) + s("BP") + struct.pack("<i", BINSIZE)
    footer += struct.pack("<q", norm_position) + struct.pack("<i", len(norm))
    parts.append(footer)

    data = bytearray(b"".join(parts))
    data[8:16] = struct.pack("<q", footer_position)
    return bytes(data)


class FakeReader:
    """A range reader over bytes in memory, counting what it was asked for."""

    def __init__(self, data: bytes) -> None:
        self.data = data
        self.requests = 0
        self.bytes = 0

    def read(self, offset: int, length: int) -> bytes:
        self.requests += 1
        out = self.data[offset : offset + length]
        self.bytes += len(out)
        return out


@pytest.fixture
def matrix():
    data = synthetic_hic()
    reader = FakeReader(data)
    return hc.HicMatrix("memory://synthetic.hic", reader=reader), data


def test_header_and_index_are_read_by_range(matrix):
    m, data = matrix
    assert m.version == 8
    assert m.chroms == CHROMS
    assert m.resolutions == [BINSIZE]
    m.read_footer(BINSIZE)
    assert m.master == {"1_1": (m.master["1_1"][0], m.master["1_1"][1])}
    assert m.norm_types(1, BINSIZE) == ["KR"]
    # the header and the footer are a handful of ranges, not a download: the synthetic file is
    # smaller than one chunk, so what is checked here is the number of requests, not their size
    assert m.reader.requests <= 6 and len(data) > 0


def test_contact_is_balanced_and_over_expected(matrix):
    m, _ = matrix
    m.read_footer(BINSIZE)
    c = m.contact("chr1", 0, 0)
    assert c["raw"] == 500 and c["norm"] == "KR" and c["same_bin"] is True
    assert c["observed"] == pytest.approx(500 / (2.0 * 2.0))
    assert c["expected"] == pytest.approx(EXPECTED[0] / SCALE_FACTOR)
    assert c["oe"] == pytest.approx(125.0 / 50.0)
    # a cell four bins away, with its own two weights and the expected value of that distance
    c = m.contact("chr1", 1_000, 4 * BINSIZE + 10)
    assert c["raw"] == 40 and c["same_bin"] is False
    assert c["observed"] == pytest.approx(40 / (2.0 * 0.5))
    assert c["oe"] == pytest.approx(40.0 / (EXPECTED[4] / SCALE_FACTOR))


def test_a_cell_with_no_record_reads_as_zero_contact(matrix):
    m, _ = matrix
    m.read_footer(BINSIZE)
    c = m.contact("chr1", 5 * BINSIZE, 8 * BINSIZE)
    assert c["raw"] == 0.0 and c["observed"] == 0.0 and c["oe"] == 0.0


def test_a_bin_the_balancing_dropped_has_no_contact(matrix):
    m, _ = matrix
    m.read_footer(BINSIZE)
    assert math.isnan(m.norm_vector("KR", 1, BINSIZE)[60])
    assert m.contact("chr1", 10 * BINSIZE, 60 * BINSIZE) is None  # KR is NaN there, no other vector
    with pytest.raises(LookupError):
        m.contact("chr7", 0, 0)


def test_the_matrix_is_symmetric_in_its_arguments(matrix):
    m, _ = matrix
    m.read_footer(BINSIZE)
    a = m.contact("chr1", 0, 4 * BINSIZE)
    b = m.contact("chr1", 4 * BINSIZE, 0)
    assert a == b


def test_block_layouts_round_trip():
    records = {(3, 5): 11, (4, 5): 12, (3, 9): 13}
    assert hc._decode_block(zlib.decompress(block_bytes(records, 0, 0)), 8) == {
        k: float(v) for k, v in records.items()
    }
    # version 9 adds two flag bytes for the width of the bin indices; here both are 32-bit
    v9 = struct.pack("<i", 2) + struct.pack("<ii", 100, 200) + bytes([1, 1, 1, 1])
    v9 += struct.pack("<i", 1) + struct.pack("<i", 3) + struct.pack("<i", 2)
    v9 += struct.pack("<i", 1) + struct.pack("<f", 1.5) + struct.pack("<i", 2) + struct.pack("<f", 2.5)
    assert hc._decode_block(v9, 9) == {(101, 203): 1.5, (102, 203): 2.5}
    # layout 2, the dense rectangle, with the -32768 hole juicer writes for an empty cell
    dense = struct.pack("<i", 3) + struct.pack("<ii", 10, 20) + bytes([0]) + bytes([2])
    dense += struct.pack("<i", 4) + struct.pack("<h", 2)
    dense += struct.pack("<hhhh", 5, -32768, 7, 8)
    assert hc._decode_block(dense, 8) == {(10, 20): 5.0, (10, 21): 7.0, (11, 21): 8.0}


def test_block_numbers_cover_both_triangles(matrix):
    m, _ = matrix
    m.read_footer(BINSIZE)
    res = m._resolution(1, BINSIZE)
    assert m._block_numbers(res, 5, 250) == [2 * BLOCK_COLUMN_COUNT + 0, 0 * BLOCK_COLUMN_COUNT + 2]


class StubMatrix:
    """Answers contacts from a dictionary, so the cache can be tested without a file."""

    def __init__(self, values):
        self.values = values
        self.reads = 0
        self.reader = FakeReader(b"")

    def contact(self, chrom, pos1, pos2, binsize=hc.RESOLUTION):
        self.reads += 1
        key = (chrom, min(pos1, pos2) // binsize, max(pos1, pos2) // binsize)
        v = self.values.get(key)
        return None if v is None else {"norm": "KR", "observed": v, "oe": v / 10, "same_bin": False}


def test_contact_source_caches_on_disk_and_asks_once(tmp_path):
    stub = StubMatrix({("chr1", 2, 40): 9.0})
    source = hc.ContactSource({"K562": "TESTACC"}, knowledge=tmp_path, readers={"K562": stub})
    assert source.contact("K562", "chr1", 10_000, 200_000)["observed"] == 9.0
    assert source.contact("K562", "chr1", 12_000, 202_000)["observed"] == 9.0  # same bins, one read
    assert stub.reads == 1
    assert source.contact("K562", "chr1", 0, 0) is None and stub.reads == 2
    assert source.contact("HepG2", "chr1", 0, 0) is None  # no matrix for this cell line
    paths = source.save()
    assert len(paths) == 1 and paths[0].name == f"K562_TESTACC_{hc.RESOLUTION}.json"
    cached = json.loads(paths[0].read_text())
    assert cached["1:2:40"]["observed"] == 9.0 and cached["1:0:0"] is None
    # a second source reads the cache and never asks the matrix
    again = hc.ContactSource({"K562": "TESTACC"}, knowledge=tmp_path, readers={"K562": stub})
    assert again.contact("K562", "chr1", 10_000, 200_000)["observed"] == 9.0
    assert stub.reads == 2 and again.misses == 0


def test_chromosome_names_normalise_both_ways():
    assert hc._chrom_key("chr1") == "1" and hc._chrom_key("1") == "1"
    assert hc._chrom_key("chrX") == "X" and hc._chrom_key("x") == "X"
