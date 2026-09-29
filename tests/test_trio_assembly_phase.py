# SPDX-License-Identifier: AGPL-3.0-or-later
"""A trio's de novo candidates placed on a phased, parent-labelled assembly of the child."""

from genomeos.genome import individuals as ind

HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\n"


def _person(tmp_path, root, name, rows):
    src = tmp_path / f"{name}.vcf"
    src.write_text(HEAD + "".join(f"chr21\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t{g}\n" for p, r, a, g in rows))
    ind.import_vcf(src, name, root=root)


def test_assembly_class_reads_each_haplotype():
    calls = {
        "alleles": {
            (100, "A", "G"): (True, False, True),
            (200, "C", "T"): (True, True, True),
            (300, "G", "A"): (False, True, False),
        },
        "positions": [100, 200, 300, 505],
    }
    assert ind.assembly_class((100, "A", "G"), calls) == "hap1"
    assert ind.assembly_class((200, "C", "T"), calls) == "both"
    assert ind.assembly_class((300, "G", "A"), calls) == "filtered"
    assert ind.assembly_class((500, "T", "C"), calls) == "nearby"  # the assembly writes 505 instead
    assert ind.assembly_class((900, "T", "C"), calls) == "absent"
    assert ind.assembly_class((900, "T", "C"), None) == "absent"


def test_read_assembly_calls_keeps_haplotype_order(tmp_path):
    vcf = tmp_path / "asm.vcf"
    vcf.write_text(
        HEAD
        + "chr21\t100\t.\tA\tG\t30\tPASS\t.\tGT\t1|0\n"
        + "chr21\t200\t.\tCT\tC,CTT\t30\tPASS\t.\tGT\t2|1\n"
        + "chr21\t300\t.\tG\tA\t30\tHET1\t.\tGT\t.|1\n"
        + "chr22\t100\t.\tA\tG\t30\tPASS\t.\tGT\t1|1\n"
    )
    c = ind.read_assembly_calls(vcf, {"chr21"})
    a = c["chr21"]["alleles"]
    assert a[(100, "A", "G")] == (True, False, True)
    assert a[(200, "CT", "C")] == (False, True, True)  # allele 1 on the second haplotype
    assert a[(200, "CT", "CTT")] == (True, False, True)  # as written: no reference, no normalising
    assert a[(300, "G", "A")] == (False, True, False)
    assert "chr22" not in c


def test_phase_trio_by_assembly_counts_candidates_and_the_control(tmp_path):
    root = tmp_path / "individuals"
    # kid: 100 from dad (control), 150 from mum (control), 300/400/500/600 in neither parent (candidates)
    _person(
        tmp_path,
        root,
        "kid",
        [
            (100, "A", "G", "0/1"),
            (150, "T", "C", "0/1"),
            (300, "G", "A", "0/1"),
            (400, "T", "C", "0/1"),
            (500, "C", "A", "0/1"),
            (600, "A", "T", "0/1"),
        ],
    )
    _person(tmp_path, root, "dad", [(100, "A", "G", "0/1")])
    _person(tmp_path, root, "mum", [(150, "T", "C", "0/1")])
    asm = tmp_path / "asm.vcf"
    asm.write_text(
        HEAD
        + "chr21\t100\t.\tA\tG\t30\tPASS\t.\tGT\t1|0\n"  # dad's allele on hap1: right
        + "chr21\t150\t.\tT\tC\t30\tPASS\t.\tGT\t1|0\n"  # mum's allele on hap1: wrong
        + "chr21\t300\t.\tG\tA\t30\tPASS\t.\tGT\t1|0\n"  # candidate on the paternal haplotype
        + "chr21\t400\t.\tT\tC\t30\tPASS\t.\tGT\t0|1\n"  # candidate on the maternal haplotype
        + "chr21\t503\t.\tG\tT\t30\tPASS\t.\tGT\t0|1\n"  # near 500 but another allele
    )
    dip = tmp_path / "dip.bed"
    dip.write_text("chr21\t0\t550\n")  # 600 lies outside the diploid assembly
    r = ind.phase_trio_by_assembly("kid", "dad", "mum", asm, dip, root=root, reference=tmp_path)
    assert r["candidates"] == {"total": 4, "outside_assembly": 1, "father": 1, "mother": 1, "nearby": 1}
    assert r["control"] == {"total": 2, "outside_assembly": 0, "right_parent": 1, "wrong_parent": 1}
    assert r["control_parent_agreement"] == 0.5 and r["candidates_paternal_share"] == 0.5
    assert r["candidates_in_benchmark_regions"] is None
