"""The proof that the blocked-gzip reference cache loses nothing.

Two halves. The first builds FASTA files here and checks the round trip
exhaustively: every byte back, thousands of random loci identical through both
readers, and a deliberately corrupted file that the same checks reject -- so the
test can fail, which is the only thing that makes the passing case worth
anything. The second reads the real reference cache when it is present and
re-proves the conversion against the hashes recorded when each chromosome was
converted: the flat `.fa` is gone, its SHA-256 is not, and the blocked file
still decompresses to exactly it.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import random
from pathlib import Path

import pytest

from genomeos.coords import Locus, Strand
from genomeos.genome import bgzf
from genomeos.genome.index import IndexedGenome, resolve_fasta, write_fai

REFERENCE = Path("data/reference")
MANIFEST = Path("data/results/reference_bgzf.json")

# Everything a real assembly can contain: soft-masked repeats in lower case, N
# runs, IUPAC ambiguity codes, a short final line, and a second record.
ALPHABET = "ACGTacgtNnRYKMSWBDHV"


def _fasta_text(seed: int, records: tuple[tuple[str, int], ...], width: int = 60) -> str:
    rng = random.Random(seed)
    out = []
    for name, length in records:
        out.append(f">{name} synthetic record for the bgzf round trip\n")
        s = "".join(rng.choice(ALPHABET) for _ in range(length))
        # A long N run and a long soft-masked run, as real chromosomes have.
        s = "N" * 500 + s[:1000] + "acgtacgt" * 250 + s[1000:]
        for i in range(0, len(s), width):
            out.append(s[i : i + width] + "\n")
    return "".join(out)


@pytest.fixture
def flat(tmp_path: Path) -> Path:
    p = tmp_path / "test.fa"
    p.write_text(_fasta_text(11, (("chrT", 40_000), ("chrU", 9_000))))
    write_fai(p)
    return p


def test_round_trip_is_byte_identical(flat: Path, tmp_path: Path) -> None:
    gz = tmp_path / "test.fa.gz"
    bgzf.compress_file(flat, gz, level=9)
    assert bgzf.is_bgzf(gz)
    ok, original, restored = bgzf.verify(flat, gz)
    assert ok, f"{original} != {restored}"
    assert bgzf.BgzfReader(gz).read() == flat.read_bytes()


def test_a_corrupted_file_is_caught(flat: Path, tmp_path: Path) -> None:
    """The proof is worth nothing if it cannot fail: change one base and it does."""
    bad = tmp_path / "bad.fa"
    raw = bytearray(flat.read_bytes())
    i = len(raw) // 2
    raw[i] = ord("A") if raw[i : i + 1] != b"A" else ord("C")
    bad.write_bytes(bytes(raw))
    gz = tmp_path / "bad.fa.gz"
    bgzf.compress_file(bad, gz, level=9)
    ok, _, _ = bgzf.verify(flat, gz)
    assert not ok


def test_still_a_gzip_file_for_every_streaming_reader(flat: Path, tmp_path: Path) -> None:
    """Genome.from_fasta and the compressor read `.fa.gz` with gzip; that must keep working."""
    gz = tmp_path / "test.fa.gz"
    bgzf.compress_file(flat, gz, level=6)
    with gzip.open(gz, "rb") as fh:
        assert fh.read() == flat.read_bytes()
    assert gzip.decompress(gz.read_bytes()) == flat.read_bytes()


def test_index_rebuilt_from_the_file_matches_the_written_one(flat: Path, tmp_path: Path) -> None:
    gz = tmp_path / "test.fa.gz"
    bgzf.compress_file(flat, gz, level=6)
    written = bgzf.gzi_path(gz).read_bytes()
    bgzf.gzi_path(gz).unlink()
    bgzf.index_file(gz)
    assert bgzf.gzi_path(gz).read_bytes() == written


def test_every_locus_reads_the_same_through_both_readers(flat: Path, tmp_path: Path) -> None:
    # The blocked copy goes in its own directory so that the flat reader is not
    # resolved away from the file it is meant to be the control for.
    gz = tmp_path / "blocked" / "test.fa.gz"
    gz.parent.mkdir()
    bgzf.compress_file(flat, gz, level=9)
    a = IndexedGenome(flat)
    b = IndexedGenome(gz)
    assert not a.blocked and b.blocked
    assert a.lengths == b.lengths
    rng = random.Random(5)
    for chrom, length in a.lengths.items():
        for _ in range(2_000):
            start = rng.randrange(length)
            end = min(length, start + rng.choice([1, 2, 7, 60, 61, 119, 200, 5_000, 20_000]))
            strand = rng.choice([Strand.PLUS, Strand.MINUS])
            locus = Locus(chrom, start, end, strand)
            assert str(a.fetch(locus)) == str(b.fetch(locus)), locus
        whole = Locus(chrom, 0, length)
        assert str(a.fetch(whole)) == str(b.fetch(whole))
    a.close()
    b.close()


def test_block_boundaries_are_not_special(flat: Path, tmp_path: Path) -> None:
    """Reads that straddle, start on and end on a block edge are the ones a bug hides in."""
    gz = tmp_path / "test.fa.gz"
    bgzf.compress_file(flat, gz, level=9)
    raw = flat.read_bytes()
    r = bgzf.BgzfReader(gz)
    edges = [0, bgzf.BLOCK_SIZE, 2 * bgzf.BLOCK_SIZE, len(raw)]
    for e in edges:
        for lo, hi in ((max(e - 10, 0), min(e + 10, len(raw))), (e, min(e + 1, len(raw)))):
            assert r.read_range(lo, hi) == raw[lo:hi], (lo, hi)
    r.close()


def test_resolve_prefers_the_blocked_file_but_keeps_the_flat_one_working(tmp_path: Path) -> None:
    flat = tmp_path / "c.fa"
    flat.write_text(_fasta_text(2, (("chrV", 3_000),)))
    assert resolve_fasta(flat) == flat
    assert resolve_fasta(tmp_path / "c.fa.gz") == flat  # asked for gz, only flat exists
    gz = tmp_path / "c.fa.gz"
    bgzf.compress_file(flat, gz, level=6)
    assert resolve_fasta(flat) == gz
    assert resolve_fasta(gz) == gz
    flat.unlink()
    assert resolve_fasta(flat) == gz


def test_a_plain_gzip_file_is_not_treated_as_blocked(tmp_path: Path) -> None:
    """A freshly downloaded UCSC `.fa.gz` is ordinary gzip: it must still be decompressed."""
    flat = tmp_path / "p.fa"
    flat.write_text(_fasta_text(3, (("chrW", 5_000),)))
    gz = tmp_path / "p.fa.gz"
    gz.write_bytes(gzip.compress(flat.read_bytes()))
    flat.unlink()
    assert not bgzf.is_bgzf(gz)
    g = IndexedGenome(gz)
    assert not g.blocked
    assert g.lengths["chrW"] > 0
    g.close()


# --- the real reference cache -------------------------------------------------


def test_a_compressed_vcf_parses_to_the_same_variants(tmp_path: Path) -> None:
    """The split HG002 files are compressed too; what matters is what iter_vcf gets back."""
    from genomeos.genome import iter_vcf

    plain = tmp_path / "S_chr21.vcf"
    rng = random.Random(9)
    lines = ["##fileformat=VCFv4.2\n", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"]
    for i in range(5_000):
        ref, alt = rng.choice(["A", "C", "G", "T", "ACGT", "AT"]), rng.choice(["A", "C", "G", "TTT"])
        gt = rng.choice(["0|1", "1|1", "0/1", "1|0"])
        lines.append(f"chr21\t{i * 7 + 1}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
    plain.write_text("".join(lines))
    gz = tmp_path / "S_chr21.vcf.gz"
    bgzf.compress_file(plain, gz, level=9)

    ok, _, _ = bgzf.verify(plain, gz)
    assert ok
    a = list(iter_vcf(plain, pass_only=False))
    b = list(iter_vcf(gz, pass_only=False))
    assert len(a) == 5_000
    assert a == b


def _manifest(key: str = "chromosomes") -> dict:
    if not MANIFEST.exists():
        return {}
    return json.loads(MANIFEST.read_text()).get(key, {})


@pytest.mark.parametrize("name", sorted(_manifest("variant_files")))
def test_converted_variant_file_still_hashes_to_the_original(name: str) -> None:
    row = _manifest("variant_files")[name]
    gz = REFERENCE / (name + ".gz")
    if not gz.exists():
        pytest.skip(f"{gz} not present in this checkout")
    assert bgzf.sha256_of_bgzf(gz) == row["sha256"], f"{name}: bytes differ from the original .vcf"


@pytest.mark.parametrize("chrom", sorted(_manifest()))
def test_converted_chromosome_still_hashes_to_the_original(chrom: str) -> None:
    """Each converted chromosome decompresses to the exact bytes of the `.fa` that was deleted."""
    row = _manifest()[chrom]
    gz = REFERENCE / f"{chrom}.fa.gz"
    if not gz.exists():
        pytest.skip(f"{gz} not present in this checkout")
    if not bgzf.is_bgzf(gz):
        # a chromosome fetched after the conversion arrives as UCSC serves it, plain gzip, and is
        # blocked when it is first read; the manifest describes what was converted, not what a fresh
        # checkout happens to hold, so this is "not converted yet" rather than "converted wrongly"
        pytest.skip(f"{gz} is plain gzip: fetched since the conversion, not yet blocked")
    h = hashlib.sha256()
    with bgzf.BgzfReader(gz) as r:
        assert len(r) == row["uncompressed_bytes"], f"{chrom}: length changed"
        pos = 0
        while pos < len(r):
            h.update(r.read_range(pos, min(pos + (1 << 24), len(r))))
            pos += 1 << 24
    assert h.hexdigest() == row["sha256"], f"{chrom}: bytes differ from the original .fa"


@pytest.mark.parametrize("chrom", sorted(_manifest())[:3])
def test_converted_chromosome_serves_loci(chrom: str) -> None:
    """Random access through the blocked file agrees with a gzip stream of the same file."""
    gz = REFERENCE / f"{chrom}.fa.gz"
    if not gz.exists():
        pytest.skip(f"{gz} not present in this checkout")
    g = IndexedGenome(gz)
    assert g.blocked
    length = g.lengths[chrom]
    rec = g.index[chrom]
    rng = random.Random(len(chrom) * 17)
    loci: list[tuple[int, int]] = []  # sorted and disjoint, so one forward pass can answer them all
    pos = 0
    while pos < length - 5_000 and len(loci) < 200:
        start = pos + rng.randrange(1, max(2, (length - 5_000) // 200))
        end = min(start + rng.randrange(1, 5_000), length)
        loci.append((start, end))
        pos = end

    def byte_span(start: int, end: int) -> tuple[int, int]:
        a = rec.offset + (start // rec.line_bases) * rec.line_bytes + start % rec.line_bases
        b = rec.offset + ((end - 1) // rec.line_bases) * rec.line_bytes + (end - 1) % rec.line_bases + 1
        return a, b

    # One forward gzip pass, which knows nothing of the .gzi, supplies the answers.
    with gzip.open(gz, "rb") as fh:
        pos = 0
        for start, end in loci:
            a, b = byte_span(start, end)
            fh.read(a - pos)
            want = fh.read(b - a)
            pos = b
            assert str(g.fetch(Locus(chrom, start, end))) == want.replace(b"\n", b"").decode().upper()
    g.close()
