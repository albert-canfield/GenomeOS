"""Telomere length from sequencing reads (task 3.2), TelSeq-style.

Reads with at least `k` TTAGGG (or CCCTAA) repeats are counted as telomeric.
Their fraction, scaled by genome size and read length, gives an estimate of
the mean telomere length per chromosome end. Inputs: FASTQ (plain or gzip) or
BAM (BGZF is gzip-compatible, so the standard library can decode it).

The method follows Ding et al. 2014 (TelSeq) as far as it can be derived: the
telomeric read fraction times the genome over its 92 ends. TelSeq's own scale
divides by the reads of comparable GC and multiplies by a constant fitted to
Southern-blot lengths; that constant is not derivable here, so the GC-corrected
ratio is reported as a dimensionless index beside the base-pair figure. Read
against GIAB's 300x BAM by ranges (docs/DATA.md).
"""

from __future__ import annotations

import gzip
import io
import json
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


TELOMERIC_GC = 0.5  # TTAGGG is three of six bases G or C, so telomeric reads sit in the 48-52% bins
GC_WINDOW = 0.02  # TelSeq's denominator: reads whose GC is within two points of a telomeric read's
GC_SAMPLE_READS = 200_000  # mapped reads read for the GC distribution (the unmapped tail is biased)


def gc_fraction(seq: str) -> float | None:
    """G+C of a read over its called bases, or None when it has none."""
    called = sum(1 for b in seq if b in "ACGTacgt")
    return (sum(1 for b in seq if b in "GCgc") / called) if called else None


def comparable_gc(seq: str, centre: float = TELOMERIC_GC, window: float = GC_WINDOW) -> bool:
    """Whether a read's GC is in the 48-52% band a telomeric read falls in (TelSeq's denominator)."""
    g = gc_fraction(seq)
    return g is not None and abs(g - centre) <= window + 1e-9  # 0.48 is inside a two-point window


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
    gc_mapped_seen = gc_mapped_hits = 0
    for chrom in chroms:  # one window inside a chromosome, away from its ends: the mapped GC distribution
        se = sequence_ends(chrom, reference)
        if se is None or chrom not in idx.by_name or gc_mapped_seen >= GC_SAMPLE_READS:
            continue
        first, last = se
        mid = (first + last) // 2
        for rd in bam.reads(chrom, mid, mid + window):
            gc_mapped_seen += 1
            gc_mapped_hits += comparable_gc(rd.seq)
            if gc_mapped_seen >= GC_SAMPLE_READS:
                break
        if progress:
            progress(
                f"GC sample {chrom}:{mid:,}: {gc_mapped_seen:,} mapped reads, {gc_mapped_hits:,} in band"
            )
    n_un = tel_un = bases_un = gc_un = 0
    for rd in bam.unmapped_sample(sample_bytes):
        n_un += 1
        bases_un += len(rd.seq)
        if is_telomeric(rd.seq, k):
            tel_un += 1
        if comparable_gc(rd.seq):
            gc_un += 1
    read_len = bases_un / n_un if n_un else 0.0
    no_coor = idx.no_coor or 0
    mapped_total = idx.mapped_total or 0
    tel_unmapped_total = (tel_un / n_un) * no_coor if n_un else 0.0
    total_reads = mapped_total + no_coor
    telomeric_total = tel_unmapped_total + mapped_tel
    fraction = telomeric_total / total_reads if total_reads else 0.0
    telomere_bp = fraction * genome_bp / ends
    # TelSeq's own denominator: reads whose GC is comparable to a telomeric read's, not every read
    gc_share_unmapped = gc_un / n_un if n_un else None
    gc_share_mapped = gc_mapped_hits / gc_mapped_seen if gc_mapped_seen else None
    # the tail is GC-biased, so the two populations are counted apart and added; when no mapped read was
    # sampled (a thin file, or a middle window with no coverage) the tail's share stands in for both
    mapped_share = gc_share_mapped if gc_share_mapped is not None else gc_share_unmapped
    gc_reads = (gc_share_unmapped or 0.0) * no_coor + (mapped_share or 0.0) * mapped_total
    gc_share = gc_reads / total_reads if total_reads else None
    # TelSeq divides by the GC-comparable reads and multiplies by a constant fitted to Southern-blot
    # lengths, which is not derivable here; the ratio itself is kept as a dimensionless index, comparable
    # between people read the same way, and the base-pair figure keeps the all-reads derivation
    gc_index = telomeric_total / gc_reads if gc_reads else None
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
        "telomeric_per_gc_comparable_read": round(gc_index, 6) if gc_index else None,
        "gc_comparable_share_of_reads": round(gc_share, 4) if gc_share is not None else None,
        "gc_comparable_share_unmapped": round(gc_share_unmapped, 4)
        if gc_share_unmapped is not None
        else None,
        "gc_comparable_share_mapped": round(gc_share_mapped, 4) if gc_share_mapped is not None else None,
        "gc_sample_mapped_reads": gc_mapped_seen,
        "gc_mapped_share_assumed_from_tail": gc_share_mapped is None,
        "gc_comparable_reads_estimated": round(gc_reads) if gc_reads else None,
        "gc_window": [round(TELOMERIC_GC - GC_WINDOW, 2), round(TELOMERIC_GC + GC_WINDOW, 2)],
        "k": k,
        "window": window,
        "bytes_fetched": bam.bytes_fetched,
        "requests": bam.requests,
        "parameter": Parameter(
            "telomere_bp_measured", round(telomere_bp), "bp", RANGE_EVIDENCE, 0.3 if n_un >= 100_000 else 0.15
        ),
        "note": "`telomere_bp` is the telomeric read fraction times the genome over its 92 ends, which is "
        "a derivation, not a calibration, and reads low against Southern blots. TelSeq instead divides by "
        "the reads of comparable GC (48 to 52%, since TTAGGG is half G or C) and multiplies by a constant "
        "fitted to blot lengths; that constant is not derivable here, so the GC correction is reported as "
        "`telomeric_per_gc_comparable_read`, a dimensionless index. Both compare people read the same way; "
        "neither is a blot length. The unmapped tail, where the pure-repeat reads sit, is sampled and "
        "scaled by the index's unmapped count, the GC shares are measured on that tail and on 200,000 "
        "mapped reads from the middle of chromosomes, and the chromosome ends are read whole",
    }


def saved_estimate(
    name: str, results: Path = Path("data/results"), individuals: Path = Path("data/individuals")
):
    """A telomere estimate already computed for this person: the range-read result under data/results (an
    open-consent person) or a telomere.json under the person's own directory. None when there is none."""
    for p in (results / f"telomere_range_{name}.json", individuals / name / "telomere.json"):
        if p.exists():
            try:
                d = json.loads(p.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            if d.get("telomere_bp") is not None:
                return {
                    "telomere_bp": float(d["telomere_bp"]),
                    "confidence": d.get("confidence"),
                    "source": str(p),
                    "note": "TelSeq-scale estimate from reads (k repeats per read), not a Southern-blot "
                    "length; calibrate the twin's attrition against it with care",
                }
    return None
