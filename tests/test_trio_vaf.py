# SPDX-License-Identifier: AGPL-3.0-or-later
"""Read allele fractions for the trio's candidates: the binomial tail, the per-call class and the
FORMAT depth reader (docs/DATA.md, "Allele fractions of the trio candidates", 2026-09-28)."""

import math

from genomeos.genome import individuals as ind


def test_binomial_tail_matches_exact_sums():
    assert math.isclose(ind.binomial_tail(0, 4), 1 / 16)
    assert math.isclose(ind.binomial_tail(2, 4), 11 / 16)
    assert math.isclose(ind.binomial_tail(4, 4, lower=False), 1 / 16)
    assert math.isclose(ind.binomial_tail(4, 4), 1.0)
    assert ind.binomial_tail(0, 0) == 1.0


def test_vaf_class_bands():
    assert ind.vaf_class(10, 10) == "shallow"
    assert ind.vaf_class(100, 100) == "half"
    assert ind.vaf_class(140, 60) == "low"  # 0.30 at depth 200
    assert ind.vaf_class(60, 140) == "high"
    assert ind.vaf_class(22, 13) == "half"  # 0.37 at depth 35: not resolvable from 0.5
    assert ind.vaf_class(125, 75) == "low"  # 0.375 at depth 200, p about 0.0002
    assert ind.vaf_class(300, 200) == "half"  # 0.40 is not below 0.40, however small p is


def test_read_allele_depths_matches_trio_keys(tmp_path):
    p = tmp_path / "d.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP:ADALL:AD\t0/1:40:18,20:9,10\n"
        "chr1\t200\t.\tC\tT,G\t50\tPASS\t.\tGT:ADALL\t1/2:1,15,17\n"
        "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT\t0/1\n"
        "chr2\t100\t.\tA\tG\t50\tPASS\t.\tGT:ADALL\t0/1:5,5\n"
    )
    d = ind.read_allele_depths(p, {"chr1"})
    assert d["chr1"][(100, "A", "G")] == (18, 20)
    assert d["chr1"][(200, "C", "G")] == (1, 17)
    assert (300, "G", "A") not in d["chr1"]
    assert "chr2" not in d
    assert ind.read_allele_depths(p, {"chr1"}, field="AD")["chr1"][(100, "A", "G")] == (9, 10)
