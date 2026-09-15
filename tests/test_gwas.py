# SPDX-License-Identifier: AGPL-3.0-or-later
"""GWAS Catalog rows parsed, the ones on elements kept, enrichment over the attribution's partitions."""

import io
import zipfile

from genomeos.attribution import gwas
from genomeos.attribution.eqtl import Intervals

HEADER = "\t".join(
    [
        "DATE ADDED TO CATALOG",
        "PUBMEDID",
        "FIRST AUTHOR",
        "DATE",
        "JOURNAL",
        "LINK",
        "STUDY",
        "DISEASE/TRAIT",
        "INITIAL SAMPLE SIZE",
        "REPLICATION SAMPLE SIZE",
        "REGION",
        "CHR_ID",
        "CHR_POS",
        "REPORTED GENE(S)",
        "MAPPED_GENE",
        "UPSTREAM_GENE_ID",
        "DOWNSTREAM_GENE_ID",
        "SNP_GENE_IDS",
        "UPSTREAM_GENE_DISTANCE",
        "DOWNSTREAM_GENE_DISTANCE",
        "STRONGEST SNP-RISK ALLELE",
        "SNPS",
        "MERGED",
        "SNP_ID_CURRENT",
        "CONTEXT",
        "INTERGENIC",
        "RISK ALLELE FREQUENCY",
        "P-VALUE",
        "PVALUE_MLOG",
        "P-VALUE (TEXT)",
        "OR or BETA",
        "95% CI (TEXT)",
        "PLATFORM [SNPS PASSING QC]",
        "CNV",
    ]
)


def _row(chrom, pos, rs, trait, gene, ctx="intron_variant", mlog="12.5", pmid="1"):
    f = [""] * 34
    f[1], f[7], f[11], f[12], f[14], f[21], f[24], f[28] = pmid, trait, chrom, pos, gene, rs, ctx, mlog
    return "\t".join(f)


TSV = (
    "\n".join(
        [
            HEADER,
            _row("21", "1000", "rs1", "Height", "A"),
            _row("21", "5000", "rs2", "Asthma", "B - C"),
            _row("21", "1;2", "rs3", "Bad", "A"),
            _row("X", "700", "rs4", "Eye colour", "D"),
            _row("", "", "rs5", "Nothing", ""),
        ]
    )
    + "\n"
)


def test_parse_rows_keeps_single_grch38_positions():
    rows = list(gwas.parse_rows(io.StringIO(TSV)))
    assert [(r[0], r[1], r[2]) for r in rows] == [
        ("chr21", 1000, "rs1"),
        ("chr21", 5000, "rs2"),
        ("chrX", 700, "rs4"),
    ]
    assert rows[1][4] == "B - C" and rows[0][6] == "12.5"


def test_distil_and_summarise(tmp_path):
    zpath = tmp_path / "gwas.zip"
    with zipfile.ZipFile(zpath, "w") as zf:
        zf.writestr("gwas-catalog-download-associations-v1.0-full.tsv", TSV)
    iv = Intervals()
    iv.add("chr21", 900, 1200, "E1")
    iv.add("chr21", 4000, 4500, "E2")  # rs2 at 5000 is outside (no margin added here)
    iv.add("chrX", 600, 800, "E3")
    iv.add("chr21", 4900, 5100, "E1" + gwas.SHIFTED)  # the control copy of E1 happens to catch rs2
    iv.freeze()
    know = tmp_path / "gwas"
    s = gwas.distil(iv, know, zpath)
    assert s["associations_read"] == 3 and s["hits"] == 3
    hits = gwas.load_hits(know)
    assert set(hits) == {"E1", "E3", "E1" + gwas.SHIFTED} and hits["E1"][0]["trait"] == "Height"
    elements = [
        {
            "key": "E1",
            "predicted": {"gene": "A", "strength": "strong"},
            "inferred": {"gene": "Z"},
            "constrained_fraction": 0.5,
        },
        {"key": "E2", "predicted": None, "inferred": {"gene": "B"}, "constrained_fraction": 0.1},
        {
            "key": "E3",
            "predicted": {"gene": "Q", "strength": "weak"},
            "inferred": {"gene": "D"},
            "constrained_fraction": 0.3,
        },
    ]
    out = gwas.summarise(elements, hits)
    assert out["with_lead_variant"] == 2 and out["fraction_with_lead_variant"] == 0.667
    assert out["shifted_control_fraction"] == 0.333 and out["enrichment_over_shifted"] == 2.0
    bp = out["by_partition"]
    assert bp["named_target"] == {"elements": 2, "with_lead_variant": 1.0}
    assert bp["no_target"] == {"elements": 1, "with_lead_variant": 0.0}
    assert bp["strong_effect"]["with_lead_variant"] == 1.0
    assert bp["constrained"] == {"elements": 2, "with_lead_variant": 1.0}
    assert (
        out["catalog_gene_is_the_predicted_target"] == 0.5
        and out["catalog_gene_is_the_inferred_target"] == 0.5
    )
    assert out["top_traits"] == {"Height": 1, "Eye colour": 1}


def test_index_adds_a_shifted_copy():
    iv = gwas.index([{"chrom": "chr1", "start": 5_000, "end": 5_200, "key": "K"}], margin=100, shift=10_000)
    assert iv.at("chr1", 4_950) == ["K"] and iv.at("chr1", 15_100) == ["K" + gwas.SHIFTED]
