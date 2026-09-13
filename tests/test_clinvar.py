# SPDX-License-Identifier: AGPL-3.0-or-later
"""ClinVar carrier screen on a distilled table and a person's file (no network)."""

import gzip

from genomeos.genome import clinvar, individuals

TABLE = """# test
#chrom\tpos\tref\talt\tgene\tsignificance\tconditions\treview\tstars\tclinvar_id\trs
chr21\t100\tA\tG\tGENE1\tPathogenic\tCondition one\treviewed_by_expert_panel\t3\t11\t123
chr21\t200\tC\tG\tGENE2\tLikely_pathogenic\tCondition two\tcriteria_provided,_single_submitter\t1\t12\t
chr21\t300\tT\tA\tGENE3\tPathogenic\tCondition three\tno_assertion_criteria_provided\t0\t13\t
"""
VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tME
chr21\t100\t.\tA\tG\t.\tPASS\t.\tGT\t1/1
chr21\t200\t.\tC\tT,G\t.\tPASS\t.\tGT\t0/2
chr21\t300\t.\tT\tA\t.\tPASS\t.\tGT\t0/0
chr21\t400\t.\tG\tA\t.\tPASS\t.\tGT\t0/1
"""


def test_screen_finds_carried_pathogenic_alleles(tmp_path):
    knowledge = tmp_path / "clinvar"
    knowledge.mkdir()
    with gzip.open(clinvar.pathogenic_path(knowledge), "wt") as fh:
        fh.write(TABLE)
    src = tmp_path / "me.vcf"
    src.write_text(VCF)
    root = tmp_path / "individuals"
    individuals.import_vcf(src, "me", root=root)
    r = clinvar.screen("me", None, root, knowledge)
    assert r["variants_scanned"] == 4 and len(r["hits"]) == 2  # 300 is 0/0, 400 is not in ClinVar
    first, second = r["hits"]
    assert first["gene"] == "GENE1" and first["zygosity"] == "homozygous" and first["stars"] == 3
    assert second["gene"] == "GENE2" and second["alt"] == "G" and second["zygosity"] == "heterozygous"
    assert r["pathogenic"] == 1 and r["likely_pathogenic"] == 1 and r["homozygous"] == 1
    assert (root / "me" / "clinvar_screen.json").exists()
