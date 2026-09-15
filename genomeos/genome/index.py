"""Random-access FASTA (task 1.2).

Writes and reads samtools-compatible `.fai` indexes so that fetch(locus) on a
250 Mb chromosome reads only the bytes it needs.

Two byte sources satisfy that index. A flat `.fa` is seeked directly. A
blocked-gzip `.fa.gz` (see `bgzf.py`) is seeked through its `.gzi`, which costs
one inflate of at most 64 KiB per block touched and saves keeping the
chromosome twice on disk. A plain, non-blocked `.gz` cannot be seeked at all,
so it is still decompressed once to a sibling `.fa` as before.

Callers name a chromosome as `chrN.fa` or `chrN.fa.gz` and do not have to know
which form is on disk: `resolve_fasta` picks whichever exists, preferring the
blocked one, and `reference_fasta` does the same from a chromosome name.
"""

from __future__ import annotations

import gzip
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import bgzf
from .sequence import Locus, Sequence, Strand


@dataclass(frozen=True, slots=True)
class FaiRecord:
    name: str
    length: int
    offset: int  # byte offset of the first base
    line_bases: int
    line_bytes: int


def decompress(path: str | Path) -> Path:
    """Return an uncompressed FASTA path, decompressing `.gz` next to it if needed."""
    path = Path(path)
    if path.suffix != ".gz":
        return path
    out = path.with_suffix("")
    if not out.exists() or out.stat().st_mtime < path.stat().st_mtime:
        with gzip.open(path, "rb") as src, open(out, "wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 20)
    return out


def resolve_fasta(path: str | Path) -> Path:
    """The file to actually read for a requested FASTA path.

    A caller may ask for `chr1.fa` when only the blocked `chr1.fa.gz` is on
    disk, or the other way round. The blocked form wins when both exist,
    because it is the one that does not need a flat copy beside it.
    """
    path = Path(path)
    if path.suffix == ".gz":
        flat, gz = path.with_suffix(""), path
    else:
        flat, gz = path, path.with_name(path.name + ".gz")
    if gz.exists() and bgzf.is_bgzf(gz):
        return gz
    if flat.exists():
        return flat
    if gz.exists():
        return gz
    return path


def reference_fasta(chrom: str, reference: str | Path = Path("data/reference")) -> Path:
    """The sequence file of one chromosome in a reference directory, in whichever form is cached.

    Returns the `.fa` path when nothing is cached, so a caller that reports a
    missing file names the form a fetch would produce.
    """
    return resolve_fasta(Path(reference) / f"{chrom}.fa")


def has_reference(chrom: str, reference: str | Path = Path("data/reference")) -> bool:
    return reference_fasta(chrom, reference).exists()


def fasta_bytes(path: str | Path) -> int:
    """How many FASTA bytes are behind a path, flat or blocked; 0 when nothing is cached.

    Callers that size a run from the file (how much memory a per-base map will
    want) need the decompressed length, not the length on disk.
    """
    p = resolve_fasta(path)
    if not p.exists():
        return 0
    if p.suffix == ".gz" and bgzf.is_bgzf(p):
        with bgzf.BgzfReader(p) as r:
            return len(r)
    return p.stat().st_size


def _open_stream(path: Path):
    """A binary stream over the FASTA bytes, whether the file is flat or gzipped."""
    return gzip.open(path, "rb") if path.suffix == ".gz" else open(path, "rb")


def write_fai(fasta: str | Path) -> Path:
    """Build `<fasta>.fai`. Requires uniform line width within each record.

    For a blocked-gzip FASTA the index is built by streaming the file and
    describes offsets into the decompressed bytes, which is exactly what
    samtools writes for a bgzipped FASTA and what the reader below seeks with.
    """
    fasta = Path(fasta)
    if not (fasta.suffix == ".gz" and bgzf.is_bgzf(fasta)):
        fasta = decompress(fasta)
    fai = fasta.with_name(fasta.name + ".fai")
    records: list[FaiRecord] = []
    name = None
    length = offset = line_bases = line_bytes = 0
    pos = 0
    with _open_stream(fasta) as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    records.append(FaiRecord(name, length, offset, line_bases, line_bytes))
                name = line[1:].split()[0].decode()
                length = 0
                offset = pos + len(line)
                line_bases = line_bytes = 0
            else:
                bases = len(line.rstrip(b"\r\n"))
                if line_bases == 0:
                    line_bases, line_bytes = bases, len(line)
                length += bases
            pos += len(line)
    if name is not None:
        records.append(FaiRecord(name, length, offset, line_bases, line_bytes))
    with open(fai, "w") as out:
        for r in records:
            out.write(f"{r.name}\t{r.length}\t{r.offset}\t{r.line_bases}\t{r.line_bytes}\n")
    return fai


def read_fai(fai: str | Path) -> dict[str, FaiRecord]:
    out: dict[str, FaiRecord] = {}
    for line in Path(fai).read_text().splitlines():
        name, length, offset, lb, lw = line.split("\t")[:5]
        out[name] = FaiRecord(name, int(length), int(offset), int(lb), int(lw))
    return out


class IndexedGenome:
    """Genome-like object serving fetch(locus) from disk via a .fai index.

    Satisfies the same contract as Genome.fetch so callers do not care which
    store they are talking to. Sequences are never loaded whole.
    """

    def __init__(self, fasta: str | Path, name: str | None = None) -> None:
        path = resolve_fasta(fasta)
        self.blocked = path.suffix == ".gz" and bgzf.is_bgzf(path)
        self.path = path if self.blocked else decompress(path)
        fai = self.path.with_name(self.path.name + ".fai")
        if not fai.exists():
            # A flat sibling's index describes the same decompressed bytes, so
            # it is reused rather than rebuilt by a second pass over the file.
            flat_fai = self.path.with_name(self.path.with_suffix("").name + ".fai")
            if self.blocked and flat_fai.exists():
                shutil.copyfile(flat_fai, fai)
            else:
                write_fai(self.path)
        self.index = read_fai(fai)
        self.name = name or self.path.name
        # Long-lived handle, closed by close(); a blocked file is seeked through
        # its .gzi instead of by the operating system.
        if self.blocked:
            self._fh: bgzf.BgzfReader | Any = bgzf.BgzfReader(self.path)
        else:
            self._fh = open(self.path, "rb")  # noqa: SIM115

    def close(self) -> None:
        self._fh.close()

    def __enter__(self) -> IndexedGenome:
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    @property
    def lengths(self) -> dict[str, int]:
        return {k: v.length for k, v in self.index.items()}

    def __len__(self) -> int:
        return sum(self.lengths.values())

    def fetch(self, locus: Locus) -> Sequence:
        rec = self.index[locus.chrom]
        if locus.end > rec.length:
            raise IndexError(f"{locus} exceeds {locus.chrom} length {rec.length}")
        if locus.length == 0:
            return Sequence("")
        start_byte = (
            rec.offset + (locus.start // rec.line_bases) * rec.line_bytes + locus.start % rec.line_bases
        )
        end_byte = (
            rec.offset
            + ((locus.end - 1) // rec.line_bases) * rec.line_bytes
            + (locus.end - 1) % rec.line_bases
            + 1
        )
        if self.blocked:
            raw = self._fh.read_range(start_byte, end_byte)
        else:
            self._fh.seek(start_byte)
            raw = self._fh.read(end_byte - start_byte)
        seq = Sequence(raw.replace(b"\n", b"").replace(b"\r", b"").decode())
        return seq.reverse_complement() if locus.strand is Strand.MINUS else seq
