# SPDX-License-Identifier: AGPL-3.0-or-later
"""A person's VCF becomes a local individual that lookup and the report can see."""

import gzip

from genomeos.genome import individuals as ind

VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001
chr21\t100\trs1\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:30
chr21\t200\t.\tC\tG\t50\t.\t.\tGT\t1|1
chr21\t300\t.\tT\tA\t10\tLowQual\t.\tGT\t0/1
22\t400\t.\tG\tT\t50\tPASS\t.\tGT\t0/1
chrUn_KI270302v1\t10\t.\tA\tC\t50\tPASS\t.\tGT\t0/1
MT\t750\t.\tA\tG\t50\tPASS\t.\tGT\t1/1
"""


def test_import_lists_and_looks_up(tmp_path, monkeypatch):
    src = tmp_path / "me.vcf.gz"
    with gzip.open(src, "wt") as fh:
        fh.write(VCF)
    root = tmp_path / "individuals"
    m = ind.import_vcf(src, "me", "test calls", root=root)
    assert m["variants"] == 4 and m["chromosomes"] == {"chr21": 2, "chr22": 1, "chrM": 1}
    assert m["skipped_filtered"] == 1 and m["skipped_other_contigs"] == 1 and m["phased"] == 1
    assert m["sample_column"] == "NA00001"
    names = [p["name"] for p in ind.list_individuals(root)]
    assert names[0] == "HG002" and "me" in names
    assert ind.vcf_path("me", "chr21", root) is not None and ind.vcf_path("me", "chr1", root) is None
    rows = list(ind.rows_in(ind.vcf_path("me", "chr21", root), 150, 250))
    assert len(rows) == 1 and rows[0][1] == "200" and rows[0][9] == "1|1"
    # carriers see the imported person (HG002 is absent for this chromosome in the test root)
    monkeypatch.setattr(ind, "ROOT", root)
    from genomeos.genome import lookup as lk

    monkeypatch.setattr(lk, "INDIVIDUALS", {})
    hits = {c["individual"]: c for c in lk.carriers("chr21", 100, "A", "G")}
    assert hits["me"]["carries"] is True and hits["me"]["genotype"] == "0/1"
    assert lk.carriers("chr21", 100, "A", "T")[0]["carries"] is False
    # names are validated, the test human is reserved, a second import needs replace
    import pytest

    with pytest.raises(ValueError):
        ind.import_vcf(src, "HG002", root=root)
    with pytest.raises(ValueError):
        ind.import_vcf(src, "bad name!", root=root)
    with pytest.raises(FileExistsError):
        ind.import_vcf(src, "me", root=root)
    assert ind.import_vcf(src, "me", root=root, replace=True)["variants"] == 4
    assert ind.remove("me", root) is True and ind.remove("me", root) is False
    assert ind.remove("HG002", root) is False


def test_normalise_chrom():
    assert ind.normalise_chrom("21") == "chr21" and ind.normalise_chrom("chrX") == "chrX"
    assert ind.normalise_chrom("MT") == "chrM" and ind.normalise_chrom("M") == "chrM"
    assert ind.normalise_chrom("chr21_random") is None and ind.normalise_chrom("HLA-A") is None


def test_verdict_reads_the_mismatch_rate():
    assert ind.verdict({"variants": 0, "reference_mismatches": 0}) == "no variants to check"
    assert ind.verdict({"variants": 10_000, "reference_mismatches": 12}).startswith("matches GRCh38")
    assert ind.verdict({"variants": 10_000, "reference_mismatches": 400}).startswith("partly disagrees")
    assert ind.verdict({"variants": 10_000, "reference_mismatches": 6_000}).startswith("does not match")


def test_dossier_reads_what_is_stored(tmp_path):
    import gzip
    import json

    src = tmp_path / "me.vcf.gz"
    with gzip.open(src, "wt") as fh:
        fh.write(VCF)
    root = tmp_path / "individuals"
    ind.import_vcf(src, "me", "test calls", root=root)
    md = ind.dossier("me", root)
    assert md.startswith("# me") and "4 PASS variants" in md and "not run" in md
    hit = {
        "chrom": "chr21",
        "pos": 100,
        "ref": "A",
        "alt": "G",
        "gene": "G1",
        "genotype": "0/1",
        "zygosity": "heterozygous",
        "significance": "Pathogenic",
        "stars": 2,
        "conditions": "C",
    }
    (root / "me" / "clinvar_screen.json").write_text(
        json.dumps(
            {
                "hits": [hit],
                "variants_scanned": 4,
                "chromosomes": ["chr21"],
                "pathogenic": 1,
                "likely_pathogenic": 0,
                "homozygous": 0,
                "two_stars_or_more": 1,
                "evidence": "curated: test",
                "note": "research",
            }
        )
    )
    md = ind.dossier("me", root)
    assert "| chr21:100 A>G | G1 | 0/1 (heterozygous) | Pathogenic | 2 | C |" in md
    assert "**Truncating variants** — not run" in md
