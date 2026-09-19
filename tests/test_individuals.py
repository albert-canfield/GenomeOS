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


def test_missense_rank_reads_uniprot_annotation():
    feats_site = [{"type": "Binding site", "description": "Zn", "start": 5, "end": 5}]
    feats_dom = [{"type": "Domain", "description": "Kinase", "start": 1, "end": 300}]
    assert ind.site_class(feats_site) == "site" and ind.site_class(feats_dom) == "domain"
    assert ind.site_class([{"type": "Chain", "description": "x", "start": 1, "end": 9}]) == "none"
    rows = [
        {"gene": "B", "pos": 1, "site": "none", "zygosity": "homozygous"},
        {"gene": "A", "pos": 2, "site": "domain", "zygosity": "heterozygous"},
        {"gene": "C", "pos": 3, "site": "site", "zygosity": "heterozygous"},
        {"gene": "C", "pos": 4, "site": "site", "zygosity": "homozygous"},
    ]
    assert [(m["gene"], m["pos"]) for m in ind.rank_missense(rows)] == [
        ("C", 4),
        ("C", 3),
        ("A", 2),
        ("B", 1),
    ]


def test_trio_counts_inheritance(tmp_path):
    root = tmp_path / "individuals"

    def person(name, rows):
        src = tmp_path / f"{name}.vcf"
        src.write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\n"
            + "".join(f"chr21\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t{g}\n" for p, r, a, g in rows)
        )
        ind.import_vcf(src, name, root=root)

    person(
        "kid",
        [(100, "A", "G", "0/1"), (200, "C", "T", "1/1"), (300, "G", "A", "0/1"), (400, "T", "C", "1/1")],
    )
    person("dad", [(100, "A", "G", "0/1"), (200, "C", "T", "0/1"), (400, "T", "C", "0/1")])
    person("mum", [(200, "C", "T", "1/1")])
    r = ind.trio("kid", "dad", "mum", root=root)
    t = r["totals"]
    # 100 from dad, 200 from both, 300 in neither (de novo candidate), 400 homozygous with mum lacking it
    assert (
        t["child_variants"],
        t["inherited"],
        t["in_both_parents"],
        t["de_novo_candidates"],
        t["mendelian_errors"],
    ) == (4, 2, 1, 1, 1)
    assert r["de_novo_candidates"][0]["pos"] == 300 and r["mendelian_errors"][0]["pos"] == 400
    assert (root / "kid" / "trio_dad_mum.json").exists()
    # trusted regions: mum was only called on 150-250, so 300 and 400 are no-calls there, not events
    bed = tmp_path / "mum.bed"
    bed.write_text("chr21\t149\t250\n")
    assert ind.import_regions("mum", bed, root) == {"chr21": 1}
    bed_dad = tmp_path / "dad.bed"
    bed_dad.write_text("chr21\t0\t1000\n")
    ind.import_regions("dad", bed_dad, root)
    r = ind.trio("kid", "dad", "mum", root=root)
    t = r["totals"]
    assert (t["de_novo_candidates"], t["mendelian_errors"], t["outside_a_parent_region"]) == (0, 0, 2)
    assert r["regions"].startswith("trusted regions")


def test_normalise_variant_makes_the_three_writings_of_one_indel_compare():
    from genomeos.genome.individuals import normalise_variant

    # TGG>TGGG, T>TG and TGG>T at one position: an insertion of G (twice) and a deletion of GG
    ref_seq = "ACGTTGGA"  # 1-based: A1 C2 G3 T4 T5 G6 G7 A8
    base_at = lambda p: ref_seq[p - 1]  # noqa: E731
    assert normalise_variant(4, "TGG", "TGGG", base_at) == normalise_variant(4, "T", "TG", base_at)
    assert normalise_variant(4, "TGG", "T", base_at) == (4, "TGG", "T")
    # a deletion inside a repeat written at a later unit left-aligns to the first
    seq = "CAGTGTGTGTA"  # C1 A2 G3 T4 G5 T6 G7 T8 G9 T10 A11
    ba = lambda p: seq[p - 1]  # noqa: E731
    assert normalise_variant(6, "TGT", "T", ba) == normalise_variant(4, "TGT", "T", ba) == (2, "AGT", "A")
    # SNVs and identical alleles are untouched
    assert normalise_variant(5, "G", "A", ba) == (5, "G", "A") and normalise_variant(5, "G", "G", ba) == (
        5,
        "G",
        "G",
    )
