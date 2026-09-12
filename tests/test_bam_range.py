# SPDX-License-Identifier: AGPL-3.0-or-later
"""An indexed BAM read by ranges: header, index pseudo-bins, BGZF blocks and records, on a synthetic file."""

import gzip
import struct
import zlib

from genomeos.genome import bam_range as br


def bgzf_block(payload: bytes) -> bytes:
    """One BGZF block: a gzip member whose extra field carries the block size."""
    comp = zlib.compressobj(6, zlib.DEFLATED, -15)
    body = comp.compress(payload) + comp.flush()
    bsize = 12 + 6 + len(body) + 8 - 1
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", bsize)
    return header + body + struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload))


def record(ref_id: int, pos: int, flag: int, seq: str) -> bytes:
    packed = bytearray()
    for i in range(0, len(seq), 2):
        a = br.SEQ.index(seq[i])
        b = br.SEQ.index(seq[i + 1]) if i + 1 < len(seq) else 0
        packed.append((a << 4) | b)
    name = b"r\x00"
    body = struct.pack("<iiBBHHHiiii", ref_id, pos, len(name), 0, 4680, 0, flag, len(seq), -1, -1, 0)
    body += name + bytes(packed) + b"\xff" * len(seq)
    return struct.pack("<i", len(body)) + body


def build(tmp_path):
    text = b"@HD\tVN:1.6\tSO:coordinate\n"
    head = b"BAM\x01" + struct.pack("<i", len(text)) + text + struct.pack("<i", 2)
    for name, length in (("chrA", 40_000), ("chrB", 20_000)):
        head += struct.pack("<i", len(name) + 1) + name.encode() + b"\x00" + struct.pack("<i", length)
    tel = "TTAGGG" * 8 + "AC"  # 50 bases, telomeric at k=7
    plain = "ACGT" * 12 + "AC"
    block_a = record(0, 100, 0, plain) + record(0, 39_500, 0, tel) + record(0, 39_600, 0, plain)
    block_b = record(1, 10, 0, plain)
    block_u = record(-1, -1, 4, tel) + record(-1, -1, 4, plain) + record(-1, -1, 4, plain)
    b0, b1, b2, b3 = bgzf_block(head), bgzf_block(block_a), bgzf_block(block_b), bgzf_block(block_u)
    eof = bgzf_block(b"")
    bam = b0 + b1 + b2 + b3 + eof
    off_a, off_b, off_u = len(b0), len(b0) + len(b1), len(b0) + len(b1) + len(b2)

    # index: chrA's linear index points every 16 kb window at block_a, chrB's at block_b;
    # the pseudo-bins carry the read counts
    def ref_index(n_intv, voff, n_mapped, chunk_end):
        bai = struct.pack("<i", 2)  # two bins: a real one and the pseudo-bin
        bai += struct.pack("<Ii", 4681, 1) + struct.pack("<QQ", voff, chunk_end)
        bai += (
            struct.pack("<Ii", br.PSEUDO_BIN, 2)
            + struct.pack("<QQ", voff, chunk_end)
            + struct.pack("<QQ", n_mapped, 0)
        )
        bai += struct.pack("<i", n_intv) + struct.pack(f"<{n_intv}Q", *([voff] * n_intv))
        return bai

    bai = b"BAI\x01" + struct.pack("<i", 2)
    bai += ref_index(3, off_a << 16, 3, off_b << 16)
    bai += ref_index(2, off_b << 16, 1, off_u << 16)
    bai += struct.pack("<Q", 3)  # n_no_coor
    p = tmp_path / "t.bam"
    p.write_bytes(bam)
    (tmp_path / "t.bam.bai").write_bytes(bai)
    return p


def test_index_and_ranges(tmp_path):
    p = build(tmp_path)
    bam = br.RemoteBam(str(p), tmp_path / "t.bam.bai")
    assert [r.name for r in bam.index.refs] == ["chrA", "chrB"]
    assert bam.index.mapped_total == 4 and bam.index.no_coor == 3
    end = list(bam.reads("chrA", 30_000, 40_000))
    assert [r.pos for r in end] == [39_500, 39_600]
    assert sum(r.seq.count("TTAGGG") >= 7 for r in end) == 1
    assert [r.pos for r in bam.reads("chrB", 0, 20_000)] == [10]
    un = list(bam.unmapped_sample(10_000))
    assert len(un) == 3 and all(r.flag & 4 for r in un) and un[0].seq.startswith("TTAGGG")


def test_decode_bgzf_stops_at_a_partial_block():
    blocks = bgzf_block(b"hello") + bgzf_block(b"world")
    data, used = br.decode_bgzf(blocks[:-5])
    assert data == b"hello" and used == len(bgzf_block(b"hello"))
    assert gzip.decompress(bgzf_block(b"x")) == b"x"


def test_estimate_remote_on_the_synthetic_bam(tmp_path, monkeypatch):
    from genomeos.genome import telomere

    p = build(tmp_path)
    monkeypatch.setattr(
        telomere, "sequence_ends", lambda chrom, reference=None: (0, 40_000) if chrom == "chrA" else None
    )
    r = telomere.estimate_remote(
        str(p), tmp_path / "t.bam.bai", ["chrA", "chrZ"], window=10_000, sample_bytes=10_000, k=7
    )
    assert r["reads_total_from_index"] == 7 and r["reads_unmapped"] == 3
    assert r["unmapped_sampled"] == 3 and r["unmapped_sample_telomeric"] == 1
    q = next(e for e in r["ends"] if e["end"] == "q")
    assert q["reads"] == 2 and q["telomeric_reads"] == 1
    # telomeric: 1 (unmapped, scaled 3/3) + 1 mapped = 2 of 7 reads
    assert r["telomeric_reads_estimated"] == 2 and abs(r["fraction"] - 2 / 7) < 1e-9
    assert r["telomere_bp"] == round(2 / 7 * telomere.GENOME_BP / telomere.CHROMOSOME_ENDS)
