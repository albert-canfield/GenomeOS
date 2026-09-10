"""Task 3.2: telomere length from reads, validated on synthetic reads with a known answer."""

import gzip
import random

from genomeos.genome.telomere import estimate, is_telomeric, iter_fastq


def _synthetic_reads(n: int, telomere_bp: float, read_len: int = 100, seed: int = 1):
    """Reads sampled uniformly from a genome whose 92 ends carry `telomere_bp` of TTAGGG."""
    rng = random.Random(seed)
    genome = 3_100_000_000
    tel_total = telomere_bp * 92
    for _ in range(n):
        if rng.random() < tel_total / genome:
            yield "TTAGGG" * (read_len // 6)
        else:
            yield "".join(rng.choice("ACGT") for _ in range(read_len))


def test_is_telomeric():
    assert is_telomeric("TTAGGG" * 10)
    assert is_telomeric("CCCTAA" * 8 + "ACGT")
    assert not is_telomeric("TTAGGG" * 3 + "ACGTACGT" * 10)


def test_estimate_recovers_known_length():
    est = estimate(_synthetic_reads(400_000, telomere_bp=8000))
    assert est.reads == 400_000
    assert abs(est.telomere_bp - 8000) / 8000 < 0.35  # sampling noise at this depth
    assert est.parameter.evidence.source.startswith("Ding")


def test_fastq_reader(tmp_path):
    fq = tmp_path / "r.fq.gz"
    with gzip.open(fq, "wt") as fh:
        for i, s in enumerate(["ACGT" * 25, "TTAGGG" * 16]):
            fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")
    reads = list(iter_fastq(fq))
    assert len(reads) == 2 and is_telomeric(reads[1])
