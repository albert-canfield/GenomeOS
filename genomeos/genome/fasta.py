"""Streaming FASTA reader. Handles plain and gzip files without loading the
whole file into memory before the first record is available."""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from pathlib import Path

from .sequence import Sequence


def _open_text(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str, Sequence]]:
    """Yield (name, description, Sequence) for each record.

    `name` is the first word of the header line; `description` is the rest.
    """
    name: str | None = None
    desc = ""
    chunks: list[str] = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, Sequence("".join(chunks))
                header = line[1:].strip()
                name, _, desc = header.partition(" ")
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, desc, Sequence("".join(chunks))


def read_fasta(path: str | Path) -> dict[str, Sequence]:
    return {name: seq for name, _, seq in iter_fasta(path)}
