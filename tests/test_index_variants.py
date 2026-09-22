"""Tasks 1.2 and 1.3: random-access FASTA and variant application."""

import random
import time
from pathlib import Path

import pytest

from genomeos.genome import (
    Genome,
    IndexedGenome,
    Locus,
    Sequence,
    Strand,
    Variant,
    apply_variants,
    iter_vcf,
    write_fai,
)

CHR21 = Path("data/reference/chr21.fa.gz")
VCF = Path("data/reference/HG002_GRCh38_benchmark.vcf.gz")


def test_fai_roundtrip_on_demo(tmp_path):
    fa = tmp_path / "d.fa"
    fa.write_text(">a desc\nACGTACGTAC\nGTACGT\n>b\nTTTT\nGG\n")
    write_fai(fa)
    g = IndexedGenome(fa)
    assert g.lengths == {"a": 16, "b": 6}
    assert str(g.fetch(Locus("a", 8, 12))) == "ACGT"
    assert str(g.fetch(Locus("a", 0, 16))) == "ACGTACGTACGTACGT"
    assert str(g.fetch(Locus("b", 3, 6, Strand.MINUS))) == "CCA"
    g.close()


@pytest.mark.skipif(not CHR21.exists(), reason="fetch chr21 first")
def _random_loci(n: int, count: int) -> list[Locus]:
    rng = random.Random(0)
    return [
        Locus("chr21", s, s + rng.randint(1, 300)) for s in (rng.randrange(0, n - 400) for _ in range(count))
    ]


def test_indexed_fetch_matches_memory():
    """What an index is for: the same sequence as reading the chromosome, at any offset.

    Split from a timing assertion on 2026-09-21. The two were one test,
    `..._matches_memory_and_is_fast`, which failed twice in one day inside the pre-push hook and
    passed alone — so a correctness check that has never been wrong reported the build broken because
    a wall clock moved. Two questions under one name fail as one.
    """
    mem = Genome.from_fasta(CHR21)
    idx = IndexedGenome(CHR21)
    n = mem.chromosomes["chr21"].length
    assert idx.lengths["chr21"] == n

    loci = _random_loci(n, 10_000)
    got = [idx.fetch(x) for x in loci]

    for x, s in zip(loci[:200], got[:200], strict=True):
        assert s == mem.fetch(x)
    idx.close()


def test_the_index_beats_a_full_parse_only_below_a_crossover():
    """The timing claim, measured rather than asserted — and the old bound was not true at its own scale.

    `assert elapsed < 1.0` for 10,000 fetches had **1.46x headroom**: 0.686s typical against a 1.0s
    bound, so any machine 46% busier failed it. That is not a flaky test, it is a bound set too tight
    to hold, and raising the number would have been moving a bar to make a failure go away.

    Measuring the alternative is what settles it. On chr21, one full parse is about 0.34s and 10,000
    indexed fetches about 0.69s, so **at ten thousand random lookups the index is slower than reading
    the whole chromosome once**. Its win is memory, and in time only below a crossover of roughly five
    thousand lookups. So the assertion is made where it is true, relatively, on the same machine in
    the same run: both arms slow together under load, which is what the wall clock could not do.
    """
    t0 = time.perf_counter()
    Genome.from_fasta(CHR21)
    full_parse = time.perf_counter() - t0

    idx = IndexedGenome(CHR21)
    loci = _random_loci(idx.lengths["chr21"], 1_000)
    t0 = time.perf_counter()
    for x in loci:
        idx.fetch(x)
    thousand = time.perf_counter() - t0
    idx.close()

    assert thousand < full_parse, (
        f"1,000 indexed fetches took {thousand:.3f}s against {full_parse:.3f}s to parse the whole "
        "chromosome: below the crossover the index must save the read, or it is buying nothing"
    )


def test_apply_variants_snv_and_indels():
    ref = Sequence("ACGTACGTAC")
    vs = [
        Variant("c", 1, "C", ("G",), gt=(1, 1)),  # hom SNV
        Variant("c", 4, "A", ("ATT",), gt=(0, 1)),  # het insertion on hap2
        Variant("c", 7, "TAC", ("T",), gt=(1, 0), phased=True),  # deletion on hap1
    ]
    h1, s1 = apply_variants(ref, vs, 0)
    h2, s2 = apply_variants(ref, vs, 1)
    assert str(h1) == "AGGTACGT" and s1["applied"] == 2 and s1["length_delta"] == -2
    assert str(h2) == "AGGTATTCGTAC" and s2["applied"] == 2 and s2["length_delta"] == 2
    # reference mismatch is skipped, not applied
    _, s3 = apply_variants(ref, [Variant("c", 0, "G", ("T",), gt=(1, 1))], 0)
    assert s3["skipped_ref_mismatch"] == 1 and s3["applied"] == 0


@pytest.mark.skipif(not (CHR21.exists() and VCF.exists()), reason="fetch chr21 and the HG002 VCF first")
def test_hg002_chr21_haplotypes():
    variants = list(iter_vcf(VCF, {"chr21"}))
    assert len(variants) > 30_000
    ref = IndexedGenome(CHR21).fetch(Locus("chr21", 0, 46_709_983))
    hap1, st1 = apply_variants(ref, variants, 0)
    assert st1["applied"] > 20_000
    assert st1["skipped_ref_mismatch"] == 0  # the benchmark VCF matches GRCh38 exactly
    assert len(hap1) == len(ref) + st1["length_delta"]
    # a homozygous SNV must be present on the haplotype
    hom = next(v for v in variants if v.is_snv and v.gt == (1, 1))
    assert str(hap1)[hom.pos] == hom.alts[0] or st1["skipped_overlap"] > 0
