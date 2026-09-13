# SPDX-License-Identifier: AGPL-3.0-or-later
"""Feature d on an injected scorer: variants collected from the elements that reach a gene, scored, summed."""

from genomeos.predict.individual_effects import collect, haplotype_sums, intervals_for, predict_gene

REG = {
    "gene": "G1",
    "promoters": [{"id": "P1", "class": "promoter", "start": 1000, "end": 1200, "distance": 5}],
    "enhancers": [
        {"id": "E_far", "class": "enhancer", "start": 5000, "end": 5300, "distance": 40000},
        {
            "id": "E_pred",
            "class": "enhancer",
            "start": 3000,
            "end": 3200,
            "distance": 2000,
            "predicted_this_gene": True,
        },
        {"id": "E_near", "class": "enhancer", "start": 2000, "end": 2100, "distance": 900},
    ],
}
VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tME
chr21\t900\t.\tA\tG\t.\tPASS\t.\tGT\t0/1
chr21\t1100\t.\tC\tT\t.\tPASS\t.\tGT\t1|0
chr21\t2050\t.\tG\tA,C\t.\tPASS\t.\tGT\t0|1
chr21\t3100\t.\tT\tTAAAA\t.\tPASS\t.\tGT\t1/1
chr21\t5100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1
"""


def scorer(chrom, pos, ref, alt):
    table = {
        1100: [("G1", "liver", -0.6), ("G1", "brain", -0.1), ("G2", "liver", 0.3)],
        2050: [("G1", "blood", 0.25)],
        3100: [("G1", "liver", -0.15)],
        5100: [("G1", "liver", 0.02)],
    }
    return table.get(pos, [])


def test_intervals_collect_score_and_sum(tmp_path):
    iv = intervals_for(REG)
    assert [e["id"] for e in iv] == ["P1", "E_pred", "E_near", "E_far"]
    vcf = tmp_path / "me_chr21.vcf"
    vcf.write_text(VCF)
    found = collect(vcf, iv)
    assert [v["pos"] for v in found] == [
        1100,
        3100,
        2050,
        5100,
    ]  # element order, 900 is outside every element
    assert found[2]["alt"] == "A"  # first alt of a multi-allelic row
    r = predict_gene("ME", "G1", "chr21", REG, scorer, vcf, max_variants=10, cache=tmp_path / "cache")
    assert r["variants_in_elements"] == 4 and r["variants_scored"] == 4 and r["variants_moving_gene"] == 3
    assert r["strongest"]["pos"] == 1100 and r["strongest"]["effect"]["tissue"] == "liver"
    assert r["strongest"]["strongest_other_gene"]["gene"] == "G2"
    h = r["haplotypes"]
    # 1100 (1|0) -0.6 on hap1; 2050 (0|1) +0.25 on hap2; 3100 (1/1) -0.15 on both; 5100 below threshold
    assert h["hap1_log2_fold_change"] == -0.75 and h["hap2_log2_fold_change"] == 0.1
    assert h["unphased_heterozygous_effects"] == []
    # the cache answers the second time without the scorer
    calls = []
    r2 = predict_gene(
        "ME", "G1", "chr21", REG, lambda *a: calls.append(a) or [], vcf, cache=tmp_path / "cache"
    )
    assert calls == [] and r2["variants_moving_gene"] == 3


def test_unphased_heterozygotes_are_kept_apart():
    rows = [
        {"genotype": "0/1", "effect": {"log2_fold_change": 0.3}},
        {"genotype": "1/1", "effect": {"log2_fold_change": -0.2}},
        {"genotype": "0/1", "effect": {"log2_fold_change": 0.05}},
    ]
    h = haplotype_sums(rows)
    assert h["hap1_log2_fold_change"] == -0.2 and h["hap2_log2_fold_change"] == -0.2
    assert h["unphased_heterozygous_effects"] == [0.3]
