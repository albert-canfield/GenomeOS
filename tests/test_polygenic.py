# SPDX-License-Identifier: AGPL-3.0-or-later
"""A polygenic score over a person's genotypes: called, assumed reference inside regions, missing outside."""

import gzip
import io
import json

from genomeos.genome import polygenic

SCORING = (
    "###PGS CATALOG SCORING FILE\n#format_version=2.0\n#pgs_id=PGS999999\n#HmPOS_build=GRCh38\n"
    "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\thm_source\thm_rsID\thm_chr\thm_pos\thm_inferOtherAllele\n"
    "21\t100\tA\tG\t0.5\tliftover\t\t21\t1000\t\n"  # called het alt (A is alt) -> dosage 1
    "21\t200\tC\tT\t-0.25\tliftover\t\t21\t2000\t\n"  # called hom alt T, effect is ref C -> dosage 0
    "21\t300\tG\tA\t1.0\tliftover\t\t21\t3000\t\n"  # no call inside region, ref base G -> dosage 2
    "21\t400\tT\tC\t2.0\tliftover\t\t21\t5000\t\n"  # outside the trusted region -> missing
    "21\t500\tA\tC\t3.0\tliftover\t\t21\t1500\t\n"  # called but neither allele matches -> mismatch
)


def test_weights_parse_harmonised_positions():
    header, rows = polygenic.weights(io.StringIO(SCORING))
    assert header["pgs_id"] == "PGS999999" and header["HmPOS_build"] == "GRCh38"
    assert rows[0] == ("chr21", 1000, "A", "G", 0.5) and len(rows) == 5


def test_score_person_sums_dosages_with_regions(tmp_path, monkeypatch):
    root = tmp_path / "people"
    d = root / "P"
    d.mkdir(parents=True)
    (d / "manifest.json").write_text(json.dumps({"name": "P", "chromosomes": ["chr21"], "source": "t"}))
    (d / "P_chr21.vcf").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP\n"
        "chr21\t1000\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
        "chr21\t1500\t.\tG\tT\t.\tPASS\t.\tGT\t0/1\n"
        "chr21\t2000\t.\tC\tT\t.\tPASS\t.\tGT\t1/1\n"
    )
    (d / "regions_chr21.bed").write_text("chr21\t0\t4000\n")
    know = tmp_path / "pgs"
    know.mkdir()
    with gzip.open(know / "PGS999999_hmPOS_GRCh38.txt.gz", "wt") as fh:
        fh.write(SCORING)
    (know / "PGS999999.json").write_text(
        json.dumps({"id": "PGS999999", "name": "T", "trait": "test", "publication": {}})
    )
    monkeypatch.setattr(polygenic, "list_individuals", lambda r: [{"name": "P", "chromosomes": ["chr21"]}])
    monkeypatch.setattr(polygenic, "_reference_base_reader", lambda chrom: ((lambda pos: "G"), None))
    r = polygenic.score_person("P", "PGS999999", root, know)
    assert (
        r["called"] == 2 and r["assumed_reference"] == 1 and r["missing"] == 1 and r["allele_mismatch"] == 1
    )
    assert r["variants_used"] == 3 and r["coverage"] == 0.6
    # 0.5*1 + (-0.25)*0 + 1.0*2 = 2.5
    assert abs(r["raw_score"] - 2.5) < 1e-9
    assert (d / "pgs_PGS999999.json").exists()
    rows = polygenic.compare(["P"], "PGS999999", root)
    assert rows[0]["raw_score"] == 2.5 and "z_within_group" not in rows[0]
