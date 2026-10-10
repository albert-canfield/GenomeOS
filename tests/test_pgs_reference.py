# SPDX-License-Identifier: AGPL-3.0-or-later
"""The polygenic reference distribution: moments, frequency matching, the person's percentile, the panel."""

import gzip
import io
import json
import math

import pytest

from genomeos.genome import pgs_reference as ref
from genomeos.genome import polygenic


def test_moments_are_the_hardy_weinberg_independence_sums():
    m = ref.moments([(0.5, 0.2, 1000), (-1.0, 0.5, 0)])
    assert m["mean"] == pytest.approx(2 * 0.2 * 0.5 + 2 * 0.5 * -1.0)
    assert m["var"] == pytest.approx(2 * 0.2 * 0.8 * 0.25 + 2 * 0.5 * 0.5 * 1.0)
    # only the variant with an allele number adds sampling error to the mean
    assert m["se_mean"] == pytest.approx(math.sqrt(4 * 0.25 * 0.2 * 0.8 / 1000))
    assert m["variants"] == 2


def test_quantile_approximations_are_close_to_tables():
    assert ref.chi2_quantile(0.975, 10) == pytest.approx(20.483, rel=0.01)
    assert ref.chi2_quantile(0.025, 10) == pytest.approx(3.247, rel=0.02)
    assert ref.chi2_quantile(0.975, 40) == pytest.approx(59.342, rel=0.005)
    assert ref.t_quantile(0.975, 15) == pytest.approx(2.1314, rel=0.005)
    assert ref.t_quantile(0.975, 22) == pytest.approx(2.0739, rel=0.005)


def test_ensembl_alleles_trim_vcf_style_indels():
    assert ref.ensembl_alleles("A", "G") == ("A", ["G"])
    assert ref.ensembl_alleles("AT", "A") == ("T", ["-"])
    assert ref.ensembl_alleles("A", "AT") == ("-", ["T"])
    assert ref.ensembl_alleles("C", "") == ("C", [])


ENTRIES = [
    {"population": "gnomADg:asj", "allele": "G", "allele_count": 3011, "frequency": 0.8672},
    {"population": "gnomADg:asj", "allele": "T", "allele_count": 461, "frequency": 0.1328},
    {"population": "gnomADg:nfe", "allele": "G", "allele_count": 900, "frequency": 1.0},
]


def test_effect_frequency_listed_unlisted_and_absent():
    p, an = ref.effect_frequency(ENTRIES, "gnomADg:asj", "T")
    assert p == pytest.approx(0.1328) and an == round(3011 / 0.8672)
    # T is not listed for nfe while G is: it sits at 1 - 1.0 = 0
    p, an = ref.effect_frequency(ENTRIES, "gnomADg:nfe", "T")
    assert p == 0.0 and an == 900
    assert ref.effect_frequency(ENTRIES, "gnomADg:afr", "T") is None


def _record(start=11091630, alleles="G/T"):
    return {
        "mappings": [
            {"seq_region_name": "19", "start": start, "allele_string": alleles, "assembly_name": "GRCh38"}
        ],
        "populations": ENTRIES,
    }


def test_match_record_checks_position_and_allele():
    row = ("chr19", 11091630, "T", "G", -0.2)
    status, freqs = ref.match_record(row, _record())
    assert status == "matched" and freqs["gnomADg:asj"][0] == pytest.approx(0.1328)
    assert "gnomADg:afr" not in freqs
    assert ref.match_record(row, _record(start=11091631))[0] == "position_mismatch"
    assert ref.match_record(row, _record(alleles="G/C"))[0] == "effect_allele_not_at_id"
    assert ref.match_record(row, None)[0] == "no_record"


SCORING = (
    "#format_version=2.0\n#pgs_id=PGS999998\n#HmPOS_build=GRCh38\n"
    "rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\tallelefrequency_effect\t"
    "hm_source\thm_rsID\thm_chr\thm_pos\thm_inferOtherAllele\n"
    "rs1\t21\t1\tA\tG\t0.5\t0.3\tENSEMBL\trs1\t21\t1000\t\n"
    "rs2\t21\t2\tC\tT\t-0.25\t\tENSEMBL\trs2\t21\t2000\t\n"
    "\t21\t3\tG\t\t1.0\tNA\tliftover\t\t21\t3000\tA\n"
)


def test_weight_table_carries_id_and_author_frequency_in_row_order():
    header, rows = polygenic.weight_table(io.StringIO(SCORING))
    assert header["pgs_id"] == "PGS999998"
    assert rows[0] == ("chr21", 1000, "A", "G", 0.5, "rs1", 0.3)
    assert rows[1][5:] == ("rs2", None)
    assert rows[2][3] == "A" and rows[2][5:] == ("", None)  # inferred other allele, no id
    _, five = polygenic.weights(io.StringIO(SCORING))
    assert five == [r[:5] for r in rows]


def _table(rows, freqs):
    return {"variants": [[*r[:5], "", None, "matched", f] for r, f in zip(rows, freqs, strict=True)]}


def _cache(know, pgs_id, table):
    know.mkdir(parents=True, exist_ok=True)
    with gzip.open(ref.af_path(pgs_id, know), "wt") as fh:
        json.dump(table, fh)


ROWS = [("chr21", 1000 + i, "A", "G", w) for i, w in enumerate([0.5, -0.3, 0.2, 0.4, 0.1])]
FREQS = [{"gnomADg:asj": [p, 3400]} for p in (0.3, 0.6, 0.1, 0.5, 0.2)]


def _person(tmp_path, ancestry):
    root = tmp_path / "people"
    (root / "P").mkdir(parents=True)
    man = {"name": "P"}
    if ancestry:
        man["ancestry"] = ancestry
    (root / "P" / "manifest.json").write_text(json.dumps(man))
    return root


def test_person_percentile_by_hand(tmp_path):
    root = _person(tmp_path, "asj")
    know = tmp_path / "pgs"
    _cache(know, "PGS1", _table(ROWS, FREQS))
    dosage = [2, 0, 1, 1, None]  # the last variant has a frequency but no genotype
    r = ref.person_reference(
        "P", "PGS1", ROWS, dosage, root=root, knowledge=know, result=tmp_path / "none.json"
    )
    assert r["status"] == "ok" and r["group"] == "gnomADg:asj"
    comp = list(zip([0.5, -0.3, 0.2, 0.4], [0.3, 0.6, 0.1, 0.5], strict=True))
    mu = sum(2 * p * w for w, p in comp)
    sd = math.sqrt(sum(2 * p * (1 - p) * w * w for w, p in comp))
    x = 0.5 * 2 + -0.3 * 0 + 0.2 * 1 + 0.4 * 1
    z = (x - mu) / sd
    assert r["z"] == pytest.approx(z, abs=1e-4)
    assert r["percentile"] == pytest.approx(100 * ref.NormalDist().cdf(z), abs=0.01)
    v = r["variants"]
    assert v["comparable"] == 4 and v["frequency_without_genotype"] == 1 and v["without_frequency"] == 0
    iv = r["interval"]
    assert iv["ld"] is None and r["ld_note"]
    # the missing variant widens the interval around the percentile, and the envelope holds every part
    lo, hi = iv["missing_variants"]
    assert lo < r["percentile"] < hi
    assert iv["combined"][0] <= min(lo, iv["af_sampling"][0])
    assert iv["combined"][1] >= max(hi, iv["af_sampling"][1])


def test_unknown_without_a_declared_or_covered_group(tmp_path):
    know = tmp_path / "pgs"
    _cache(know, "PGS1", _table(ROWS, FREQS))
    dosage = [2, 0, 1, 1, 1]
    r = ref.person_reference("P", "PGS1", ROWS, dosage, root=_person(tmp_path, None), knowledge=know)
    assert r["status"] == "UNKNOWN" and r["percentile"] is None and "no group is borrowed" in r["reason"]
    r = ref.person_reference("P", "PGS1", ROWS, dosage, root=_person(tmp_path / "b", "nfe"), knowledge=know)
    # nfe is a gnomAD group, but this table has no nfe frequency: nothing comparable, so no variance
    assert r["status"] == "UNKNOWN" and r["variants"]["without_frequency"] == 5
    r = ref.person_reference(
        "P", "PGS1", ROWS, dosage, root=_person(tmp_path / "c", "martian"), knowledge=know
    )
    assert r["status"] == "UNKNOWN" and "not a group the source covers" in r["reason"]


def _ld(res, ratio, ci, reading):
    pooled = {"ratio": ratio, "ci": ci, "reading": reading}
    res.write_text(json.dumps({"scores": {"PGS1": {"panel": {"ld": {"pooled": pooled}}}}}))
    return res


def test_coverage_floor_and_ld_rule(tmp_path):
    root = _person(tmp_path, "asj")
    know = tmp_path / "pgs"
    _cache(know, "PGS1", _table(ROWS, FREQS))
    sparse = [2, None, None, None, None]  # 0.5 of 1.5 total |weight| comparable
    r = ref.person_reference("P", "PGS1", ROWS, sparse, root=root, knowledge=know, result=tmp_path / "x.json")
    assert r["status"] == "UNKNOWN" and "below the floor" in r["reason"]
    assert r["diagnostic"]["percentile_under_independence"] is not None
    full = [2, 0, 1, 1, 0]
    res = _ld(tmp_path / "res.json", 2.4, [1.6, 3.9], "fails")
    r = ref.person_reference("P", "PGS1", ROWS, full, root=root, knowledge=know, result=res)
    assert r["status"] == "UNKNOWN" and r["reason"].startswith("R1")
    res = _ld(tmp_path / "res.json", 1.1, [0.7, 1.8], "holds")
    r = ref.person_reference("P", "PGS1", ROWS, full, root=root, knowledge=know, result=res)
    assert r["status"] == "ok" and r["ld"]["reading"] == "holds"
    lo, hi = r["interval"]["ld"]
    assert lo <= r["percentile"] <= hi


def test_misaligned_table_is_refused(tmp_path):
    root = _person(tmp_path, "asj")
    know = tmp_path / "pgs"
    _cache(know, "PGS1", _table(ROWS[::-1], FREQS))
    r = ref.person_reference("P", "PGS1", ROWS, [1] * 5, root=root, knowledge=know)
    assert r["status"] == "UNKNOWN" and "row for row" in r["reason"]


def test_giab_trio_declares_ashkenazi(tmp_path):
    a = ref.declared_ancestry("HG002", tmp_path)
    assert a["group"] == "gnomADg:asj" and "declared" in a["source"]
    assert ref.declared_ancestry("someone", tmp_path)["group"] is None


PEOPLE = {
    f"GCA_{k:03d}": {
        "sample": f"S{k // 2}",
        "superpopulation": "AFR" if k < 24 else "AMR",
        "subpopulation": "",
    }
    for k in range(48)
}


def _genotypes(pattern):
    """(row index, {sample: dosage}) for 24 people from per-row dosage lists."""
    for i, ds in enumerate(pattern):
        yield i, {f"S{j}": d for j, d in enumerate(ds)}


def test_panel_check_ratio_matches_hand_computation_for_one_variant():
    ds = [0, 1, 1, 2, 0, 1] * 4  # 24 people, 12 per group
    rows = [("chr1", 10, "A", "G", 0.7)]
    out = ref.panel_check(rows, _genotypes([ds]), PEOPLE)
    for g, sl in (("AFR", ds[:12]), ("AMR", ds[12:])):
        n = len(sl)
        p = sum(sl) / (2 * n)
        analytic = 2 * p * (1 - p) * 0.49 * (2 * n) / (2 * n - 1)
        xs = [0.7 * d for d in sl]
        emp = sum((x - sum(xs) / n) ** 2 for x in xs) / (n - 1)
        assert out["groups"][g]["ratio"] == pytest.approx(emp / analytic, abs=1e-4)
    assert out["pooled"]["df"] == 22 and set(out["pooled"]["groups"]) == {"AFR", "AMR"}


def test_panel_check_sees_perfect_ld_as_a_multiplied_variance():
    base = [0, 1, 1, 2, 0, 1, 2, 1, 0, 1, 1, 0] * 2
    rows = [("chr1", 10 + i, "A", "G", 1.0) for i in range(8)]
    out = ref.panel_check(rows, _genotypes([base] * 8), PEOPLE)  # eight copies of one variant
    assert out["pooled"]["ratio"] > 6 and out["pooled"]["reading"] == "fails"
    # a person missing a haplotype drops the row for that group only
    gaps = [list(base) for _ in range(8)]
    gaps[0][0] = None
    out = ref.panel_check(rows, _genotypes(gaps), PEOPLE)
    assert out["groups"]["AFR"]["snv_rows_placed_in_every_haplotype"] == 7
    assert out["groups"]["AMR"]["snv_rows_placed_in_every_haplotype"] == 8


def test_panel_dosages_read_sites_blocks_and_absence(tmp_path, monkeypatch):
    d = tmp_path / "panel" / "chr1"
    d.mkdir(parents=True)
    (d / "meta.json").write_text(
        json.dumps({"assemblies": ["GCA_1.1", "GCA_2.1", "GCA_3.1", "GCA_4.1", "hs1"]})
    )
    with gzip.open(d / "sites.tsv.gz", "wt") as fh:
        fh.write("99\tG\tA=AA=\n")  # 0-based 99 = 1-based 100; X carries A/G, Y carries A/A
        fh.write("199\tC\t.=T=.\n")  # X's first haplotype absent at 200
    with gzip.open(d / "blocks.tsv.gz", "wt") as fh:
        fh.write("0\t1000\tsss.s\t0\n")  # Y's second haplotype absent across the block
    people = {
        "GCA_1": {"sample": "X", "superpopulation": "AFR", "subpopulation": ""},
        "GCA_2": {"sample": "X", "superpopulation": "AFR", "subpopulation": ""},
        "GCA_3": {"sample": "Y", "superpopulation": "AFR", "subpopulation": ""},
        "GCA_4": {"sample": "Y", "superpopulation": "AFR", "subpopulation": ""},
    }
    import genomeos.genome.individuals as ind

    monkeypatch.setattr(ind, "_reference_base_reader", lambda chrom: ((lambda pos1: "T"), None))
    rows = [
        ("chr1", 100, "A", "G", 1.0),
        ("chr1", 200, "T", "C", 1.0),
        ("chr1", 500, "T", "C", 1.0),  # no site: hg38's base where the assembly is present
        ("chr1", 600, "AT", "A", 1.0),  # not an SNV: skipped
    ]
    got = dict(ref.panel_dosages(rows, panel=tmp_path / "panel", people=people))
    assert got[0] == {"X": 1, "Y": 2}
    assert got[1] == {"X": None, "Y": 1}
    assert got[2] == {"X": 2, "Y": None}
    assert 3 not in got


def test_score_person_reports_unknown_reference_without_ancestry(tmp_path, monkeypatch):
    root = tmp_path / "people"
    d = root / "P"
    d.mkdir(parents=True)
    (d / "manifest.json").write_text(json.dumps({"name": "P", "chromosomes": ["chr21"]}))
    (d / "P_chr21.vcf").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP\n"
        "chr21\t1000\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
    )
    know = tmp_path / "pgs"
    know.mkdir()
    with gzip.open(know / "PGS999998_hmPOS_GRCh38.txt.gz", "wt") as fh:
        fh.write(SCORING)
    (know / "PGS999998.json").write_text(json.dumps({"id": "PGS999998", "publication": {}}))
    monkeypatch.setattr(polygenic, "list_individuals", lambda r: [{"name": "P", "chromosomes": ["chr21"]}])
    monkeypatch.setattr(polygenic, "_reference_base_reader", lambda chrom: ((lambda pos: "G"), None))
    r = polygenic.score_person("P", "PGS999998", root, know)
    assert r["reference"]["status"] == "UNKNOWN" and r["reference"]["percentile"] is None
    assert "not a risk" in r["reference"]["note"]


def test_describe_and_cell_word_the_percentile_and_unknown():
    ok = {
        "status": "ok",
        "percentile": 62.14,
        "group": "gnomADg:asj",
        "group_label": "Ashkenazi Jewish",
        "interval": {"combined": [48.0, 74.9]},
        "variants": {"comparable": 280, "in_score": 313, "weight_share_comparable": 0.9},
        "ld": {"reading": "holds"},
    }
    line = ref.describe(ok)
    assert "62.1" in line and "48.0 to 74.9" in line and "280 of 313" in line and "not a risk" in line
    assert ref.cell(ok) == "62.1 (48.0-74.9), asj"
    gone = {"status": "UNKNOWN", "reason": "no ancestry declared for this person, and no group is borrowed"}
    assert ref.describe(gone).startswith("percentile UNKNOWN: no ancestry")
    assert ref.cell(gone) == "UNKNOWN" and ref.cell(None) == "UNKNOWN"


def test_dossier_carries_the_percentile_column(tmp_path, monkeypatch):
    from genomeos.genome import individuals as ind

    root = tmp_path / "people"
    d = root / "P"
    d.mkdir(parents=True)
    monkeypatch.setattr(
        ind,
        "list_individuals",
        lambda r=None: [{"name": "P", "builtin": False, "chromosomes": [], "evidence": "t"}],
    )
    reference = {
        "status": "ok",
        "percentile": 10.0,
        "group": "gnomADg:nfe",
        "interval": {"combined": [5.0, 18.0]},
    }
    (d / "pgs_PGS1.json").write_text(
        json.dumps(
            {
                "score": {"id": "PGS1", "trait": "t"},
                "variants_used": 3,
                "variants_in_score": 4,
                "raw_score": 1.0,
                "reference": reference,
            }
        )
    )
    (d / "pgs_PGS2.json").write_text(
        json.dumps(
            {
                "score": {"id": "PGS2", "trait": "u"},
                "variants_used": 3,
                "variants_in_score": 4,
                "raw_score": 1.0,
            }
        )
    )
    page = ind.dossier("P", root)
    assert "| PGS1 | t |  | 3 of 4 | +1.0000 | 10.0 (5.0-18.0), nfe |" in page
    assert "| PGS2 | u |  | 3 of 4 | +1.0000 | UNKNOWN |" in page
    assert "not a risk" in page
