"""Telomere length from sequencing reads (task 3.2), TelSeq-style.

Reads with at least `k` TTAGGG (or CCCTAA) repeats are counted as telomeric.
Their fraction, scaled by genome size and read length, gives an estimate of
the mean telomere length per chromosome end. Inputs: FASTQ (plain or gzip) or
BAM (BGZF is gzip-compatible, so the standard library can decode it).

The method follows Ding et al. 2014 (TelSeq): length ≈ (telomeric reads /
reads with comparable GC) × genome_length / (46 chromosome ends), but we
normalise by *all* reads rather than a GC bin, which is simpler and adequate
for a first estimate. Validation on a real 30x BAM is pending (docs/PROGRESS.md).
"""

from __future__ import annotations

import gzip
import io
import struct
import urllib.request
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from genomeos.ir import Evidence, EvidenceKind, Parameter

RANGE_EVIDENCE = Evidence(
    EvidenceKind.INFERRED,
    "TelSeq (Ding et al. 2014) over an indexed BAM read by ranges: the unmapped tail sampled, "
    "the mapped ends read whole, totals from the index",
)
END_WINDOW = 10_000

REPEAT = "TTAGGG"
REPEAT_RC = "CCCTAA"
TELSEQ_EVIDENCE = Evidence(EvidenceKind.EXPERIMENTAL, "Ding et al. 2014, Nucleic Acids Res 42:e75 (TelSeq)")
GENOME_BP = 3_100_000_000
CHROMOSOME_ENDS = 46 * 2


@dataclass(slots=True)
class TelomereEstimate:
    reads: int
    telomeric_reads: int
    read_length: float
    telomere_bp: float
    parameter: Parameter

    @property
    def fraction(self) -> float:
        return self.telomeric_reads / self.reads if self.reads else 0.0


def is_telomeric(seq: str, k: int = 7) -> bool:
    return seq.count(REPEAT) >= k or seq.count(REPEAT_RC) >= k


def iter_fastq(path: str | Path) -> Iterator[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().strip()
            fh.readline()
            fh.readline()
            yield seq


def iter_fastq_url(url: str, chunk: int = 1 << 20) -> Iterator[str]:
    """Stream a (gzipped) FASTQ over HTTP: nothing is written to disk.

    Reads arrive in sequencer order, which is effectively random with respect
    to the genome, so the first N reads are an unbiased sample. Stop early by
    breaking out of the iterator; the connection is closed on exit."""
    req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
    with urllib.request.urlopen(req, timeout=120) as resp:  # noqa: S310  (https URL from the caller)
        raw = resp if not url.endswith(".gz") else gzip.GzipFile(fileobj=io.BufferedReader(resp, chunk))
        text = io.TextIOWrapper(raw, encoding="ascii", errors="replace")
        while True:
            header = text.readline()
            if not header:
                return
            seq = text.readline().strip()
            text.readline()
            text.readline()
            yield seq


_BAM_SEQ = "=ACMGRSVTWYHKDBN"


def iter_bam(path: str | Path) -> Iterator[str]:
    """Yield read sequences from a BAM file using only the standard library."""
    with gzip.open(path, "rb") as fh:
        magic = fh.read(4)
        if magic != b"BAM\x01":
            raise ValueError("not a BAM file")
        (l_text,) = struct.unpack("<i", fh.read(4))
        fh.read(l_text)
        (n_ref,) = struct.unpack("<i", fh.read(4))
        for _ in range(n_ref):
            (l_name,) = struct.unpack("<i", fh.read(4))
            fh.read(l_name + 4)
        while True:
            head = fh.read(4)
            if len(head) < 4:
                return
            (block_size,) = struct.unpack("<i", head)
            block = fh.read(block_size)
            l_read_name = block[8]
            n_cigar_op = struct.unpack("<H", block[12:14])[0]
            l_seq = struct.unpack("<i", block[16:20])[0]
            off = 32 + l_read_name + 4 * n_cigar_op
            packed = block[off : off + (l_seq + 1) // 2]
            chars = []
            for byte in packed:
                chars.append(_BAM_SEQ[byte >> 4])
                chars.append(_BAM_SEQ[byte & 0xF])
            yield "".join(chars[:l_seq])


def estimate(
    reads: Iterator[str],
    k: int = 7,
    genome_bp: int = GENOME_BP,
    ends: int = CHROMOSOME_ENDS,
    max_reads: int | None = None,
) -> TelomereEstimate:
    n = tel = 0
    total_len = 0
    for seq in reads:
        n += 1
        total_len += len(seq)
        if is_telomeric(seq, k):
            tel += 1
        if max_reads and n >= max_reads:
            break
    read_len = total_len / n if n else 0.0
    # telomeric bases in the genome ≈ fraction × genome; per end divide by 92
    telomere_bp = (tel / n) * genome_bp / ends if n else 0.0
    return TelomereEstimate(
        n,
        tel,
        read_len,
        telomere_bp,
        Parameter(
            "telomere_bp_measured", round(telomere_bp), "bp", TELSEQ_EVIDENCE, 0.5 if n >= 100_000 else 0.2
        ),
    )


def estimate_file(path: str | Path, **kw) -> TelomereEstimate:
    """FASTQ/BAM on disk, or an http(s) URL to a FASTQ streamed without saving it."""
    p = str(path)
    if p.startswith(("http://", "https://")):
        reads = iter_fastq_url(p)
    else:
        reads = iter_bam(p) if p.endswith(".bam") else iter_fastq(p)
    return estimate(reads, **kw)


def sequence_ends(chrom: str, reference: Path = Path("data/reference")) -> tuple[int, int] | None:
    """(first, last) non-N positions of the chromosome: where its assembled sequence starts and ends."""
    from genomeos.coords import Locus
    from genomeos.genome import IndexedGenome

    fa = reference / f"{chrom}.fa"
    if not fa.exists():
        return None
    g = IndexedGenome(str(fa))
    try:
        length = g.lengths[chrom]
        step = 100_000
        first = last = None
        for lo in range(0, length, step):
            s = str(g.fetch(Locus(chrom, lo, min(length, lo + step)))).upper()
            stripped = s.lstrip("N")
            if stripped:
                first = lo + (len(s) - len(stripped))
                break
        for hi in range(length, 0, -step):
            s = str(g.fetch(Locus(chrom, max(0, hi - step), hi))).upper()
            stripped = s.rstrip("N")
            if stripped:
                last = hi - (len(s) - len(stripped))
                break
    finally:
        g.close()
    return (first, last) if first is not None and last is not None else None


def estimate_remote(
    url: str,
    index: str | Path,
    chroms: list[str],
    reference: Path = Path("data/reference"),
    window: int = END_WINDOW,
    sample_bytes: int = 100_000_000,
    k: int = 7,
    genome_bp: int = GENOME_BP,
    ends: int = CHROMOSOME_ENDS,
    progress=None,
) -> dict:
    """TelSeq over a 600 GB BAM without downloading it: the unmapped tail sampled by ranges (where the
    pure-repeat reads are), each chromosome end read whole (where the boundary reads are), and the read
    totals taken from the index's pseudo-bins. Everything is counted, nothing kept."""
    from genomeos.genome.bam_range import RemoteBam

    bam = RemoteBam(url, index)
    idx = bam.index
    end_rows = []
    mapped_tel = 0
    for chrom in chroms:
        se = sequence_ends(chrom, reference)
        if se is None or chrom not in idx.by_name:
            continue
        first, last = se
        for side, (a, b) in (("p", (first, first + window)), ("q", (last - window, last))):
            n = tel = bases = 0
            for rd in bam.reads(chrom, a, b):
                n += 1
                bases += len(rd.seq)
                if is_telomeric(rd.seq, k):
                    tel += 1
            mapped_tel += tel
            end_rows.append(
                {
                    "chrom": chrom,
                    "end": side,
                    "window": [a, b],
                    "reads": n,
                    "telomeric_reads": tel,
                    "coverage": round(bases / window, 1) if window else None,
                    "telomeric_fraction": round(tel / n, 5) if n else None,
                }
            )
            if progress:
                progress(
                    f"{chrom} {side} end: {n:,} reads, {tel:,} telomeric; {bam.bytes_fetched / 1e6:.0f} MB"
                )
    n_un = tel_un = bases_un = 0
    for rd in bam.unmapped_sample(sample_bytes):
        n_un += 1
        bases_un += len(rd.seq)
        if is_telomeric(rd.seq, k):
            tel_un += 1
    read_len = bases_un / n_un if n_un else 0.0
    no_coor = idx.no_coor or 0
    mapped_total = idx.mapped_total or 0
    tel_unmapped_total = (tel_un / n_un) * no_coor if n_un else 0.0
    total_reads = mapped_total + no_coor
    fraction = (tel_unmapped_total + mapped_tel) / total_reads if total_reads else 0.0
    telomere_bp = fraction * genome_bp / ends
    covs = [r["coverage"] for r in end_rows if r["reads"]]
    return {
        "source": url,
        "index": str(index),
        "reads_total_from_index": total_reads,
        "reads_mapped": mapped_total,
        "reads_unmapped": no_coor,
        "unmapped_sampled": n_un,
        "unmapped_sample_telomeric": tel_un,
        "unmapped_telomeric_fraction": round(tel_un / n_un, 6) if n_un else None,
        "telomeric_reads_estimated": round(tel_unmapped_total + mapped_tel),
        "mean_read_length": round(read_len, 1),
        "mean_end_coverage": round(sum(covs) / len(covs), 1) if covs else None,
        "ends": end_rows,
        "telomere_bp": round(telomere_bp),
        "fraction": fraction,
        "k": k,
        "window": window,
        "bytes_fetched": bam.bytes_fetched,
        "requests": bam.requests,
        "parameter": Parameter(
            "telomere_bp_measured", round(telomere_bp), "bp", RANGE_EVIDENCE, 0.3 if n_un >= 100_000 else 0.15
        ),
        "note": "TelSeq counts reads with at least k TTAGGG repeats among all reads; here the unmapped tail, "
        "where the pure-repeat reads sit, is sampled and scaled by the index's unmapped count, the mapped "
        "chromosome ends are read whole, and the mapped total comes from the index; no GC normalisation, "
        "so the number is a first estimate to compare between people read the same way, not a length "
        "to quote",
    }
