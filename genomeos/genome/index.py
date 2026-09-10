"""Random-access FASTA (task 1.2).

Writes and reads samtools-compatible `.fai` indexes so that fetch(locus) on a
250 Mb chromosome reads only the bytes it needs. Compressed input is
decompressed once to a sibling `.fa` because gzip cannot be seeked.
"""

from __future__ import annotations

import gzip
import shutil
from dataclasses import dataclass
from pathlib import Path

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


def write_fai(fasta: str | Path) -> Path:
    """Build `<fasta>.fai`. Requires uniform line width within each record."""
    fasta = decompress(fasta)
    fai = fasta.with_name(fasta.name + ".fai")
    records: list[FaiRecord] = []
    with open(fasta, "rb") as fh:
        name = None
        length = offset = line_bases = line_bytes = 0
        pos = 0
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
        self.path = decompress(fasta)
        fai = self.path.with_name(self.path.name + ".fai")
        if not fai.exists():
            write_fai(self.path)
        self.index = read_fai(fai)
        self.name = name or self.path.name
        self._fh = open(self.path, "rb")  # noqa: SIM115  (long-lived handle, closed by close())

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
        self._fh.seek(start_byte)
        raw = self._fh.read(end_byte - start_byte)
        seq = Sequence(raw.replace(b"\n", b"").replace(b"\r", b"").decode())
        return seq.reverse_complement() if locus.strand is Strand.MINUS else seq
