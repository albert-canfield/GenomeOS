# SPDX-License-Identifier: AGPL-3.0-or-later
"""The demo case above data level 1: GIAB HG008 with its matched normal genome.

docs/THERAPEUTICS.md, "A demo case above data level 1", decided what level 2 requires here (the
normal's own calls, held and subtracted by the pipeline) before the case was built. These tests pin:

- the level the committed inputs reach (2, and no higher: the chain breaks at tumour RNA);
- that the report at that level refuses what needs level 3 or above, and no longer says what is only
  true without a normal;
- that the two sentences which used to assume level 1 still read word for word without a normal, so
  every level-1 output is unchanged;
- that the population rule's "likely inherited" verdict gives way to the normal, the measurement.
"""

from __future__ import annotations

import json
from pathlib import Path

from genomeos.cancer.tumour import GERMLINE_AF, grade
from genomeos.genome.variants import iter_vcf
from genomeos.therapeutics import pipeline as pl
from genomeos.therapeutics.pipeline import GERMLINE_VERDICT, PatientProfile, origin_of
from genomeos.therapeutics.report import global_limitations

ROOT = Path(__file__).resolve().parent.parent
DEMO = ROOT / "data/demo/hg008_pdac"
TUMOUR = DEMO / "HG008-T.vcf"
NORMAL = DEMO / "HG008-N-D.vcf"
PROVENANCE = DEMO / "provenance.json"
COMPARISON = DEMO / "level_comparison.json"
REPORT = DEMO / "therapeutics/therapeutic_report.txt"
ANALYSIS = DEMO / "therapeutics/therapeutic_analysis.json"

LEVEL1_NORMAL_SENTENCE = (
    "No matched normal sample means somatic status is estimated from population allele frequency; "
    "rare inherited variants cannot be separated from somatic ones."
)
LEVEL1_BASIS = "no matched normal supplied; population frequency used to set inherited variants aside"


# --- the committed inputs -------------------------------------------------------------


def test_the_demo_pair_reaches_level_two_and_no_higher():
    levels = PatientProfile("HG008-T", str(TUMOUR), normal_vcf=str(NORMAL)).levels()
    assert levels["level_reached"] == 2
    assert levels["highest_input_present"] == 2
    assert levels["inputs_present"] == ["tumour_vcf", "matched_normal"]
    # the chain breaks at tumour RNA, which HG008 publishes only as raw reads from other passages
    first_missing = next(a for a in levels["available"] if not a["present"])
    assert first_missing["input"] == "tumour_rna"


def test_the_same_tumour_file_alone_is_level_one():
    assert PatientProfile("HG008-T", str(TUMOUR)).levels()["level_reached"] == 1


def test_the_inputs_are_the_registered_sources_cut_to_one_region_set():
    p = json.loads(PROVENANCE.read_text())
    # the md5s written into docs/THERAPEUTICS.md before the build
    assert p["tumour"]["md5"] == "72e85c9ddd6c222137e4136f6a76d21f"
    assert p["normal"]["md5"] == "891836dd9da087bceea43ee3e436f1dc"
    assert p["tumour"]["url"].endswith("HG008-T--HG008-N.snv.indel.high_confidence.v7.annotated.vcf.gz")
    assert p["normal"]["url"].endswith("HG008-N.haplotypecaller.gatk.annotated.unfiltered.vcf.gz")
    assert "padded by 2 bp" in p["region_set"]["rule"]
    # what the files hold is what the provenance says was written
    assert sum(1 for _ in iter_vcf(TUMOUR, None, pass_only=False)) == p["tumour"]["records_written"] == 107
    assert sum(1 for _ in iter_vcf(NORMAL, None, pass_only=False)) == p["normal"]["records_written"] == 33714
    # every normal record carries an alternate allele, because the subtraction does not read genotypes
    assert all(any(a > 0 for a in v.gt) for v in iter_vcf(NORMAL, None, pass_only=False))
    # small enough to ship, and every file names its source and what was changed
    assert sum(f.stat().st_size for f in DEMO.rglob("*") if f.is_file()) < 20_000_000
    for f in (TUMOUR, NORMAL):
        head = f.read_text()[:2000]
        assert "National Institute of Standards and Technology" in head
        assert "Changed by GenomeOS" in head


def test_the_subtraction_count_is_the_one_the_comparison_reports():
    from genomeos.cancer.compare import somatic

    removed = 107 - len(somatic(str(NORMAL), str(TUMOUR)))
    c = json.loads(COMPARISON.read_text())
    assert c["level_2"]["removed_by_the_normal"] == removed
    assert c["level_1_same_tumour_file"]["removed_by_the_normal"] == 0


# --- the committed run ----------------------------------------------------------------


def test_the_committed_run_states_level_two_and_refuses_what_needs_more():
    c = json.loads(COMPARISON.read_text())
    two, one = c["level_2"], c["level_1_same_tumour_file"]
    assert two["level_reached"] == 2 and one["level_reached"] == 1
    assert "matched normal genome" not in two["missing_data"]
    assert "matched normal genome" in one["missing_data"]
    assert "tumour RNA-seq" in two["missing_data"]
    # at level 2 nothing is set aside by population frequency: the normal decided
    assert two["likely_germline"] == 0
    assert all(s == "somatic_candidate" for cand in two["candidates"] for s in cand["somatic_status"])
    # no candidate's expression in this tumour is claimed: that needs level 3
    assert two["candidates"] and not any(cand["expression_known"] for cand in two["candidates"])

    report = REPORT.read_text()
    assert "Data level reached: 2 of 9" in report
    assert "This analysis ran at data level 2 of 9. Conclusions that need a higher" in report
    assert "No tumour RNA or protein data was supplied" in report
    assert "A matched normal genome was supplied" in report
    assert "No matched normal sample" not in report
    analysis = json.loads(ANALYSIS.read_text())
    assert analysis["data_level"] == 2
    origins = [o for cand in analysis["therapeutic_candidates"] for o in cand["origin_variants"]]
    assert origins and all(o["somatic_basis"].startswith("absent from the matched normal") for o in origins)


def test_the_one_record_the_population_rule_sets_aside_is_somatic_by_the_normal():
    """What level 2 changed on this tumour, measured: one record, and the normal's reads agree.

    chrX:24066124 G>GT (EIF2S3, intronic) carries a gnomAD frequency of 0.8185 in its VEP annotation,
    so level 1 sets it aside as likely inherited. The normal's calls do not carry it, so level 2 keeps
    it. NYGC's own pileup of the normal reads 105 reference and 1 alternate read there.
    """
    c = json.loads(COMPARISON.read_text())
    assert c["level_1_same_tumour_file"]["likely_germline"] == 1
    assert c["level_2"]["likely_germline"] == 0
    assert c["level_2"]["removed_by_the_normal"] == 0
    tumour = {(v.chrom, v.pos + 1, v.ref, v.alts) for v in iter_vcf(TUMOUR, None, pass_only=False)}
    normal = {(v.chrom, v.pos + 1) for v in iter_vcf(NORMAL, None, pass_only=False)}
    assert ("chrX", 24066124, "G", ("GT",)) in tumour
    assert ("chrX", 24066124) not in normal


# --- the wording that used to assume level 1 -------------------------------------------


def _analysis(inputs: list[str]) -> dict:
    level = 2 if "matched_normal" in inputs else 1
    return {"data_level": {"level_reached": level, "inputs_present": inputs}, "candidates": []}


def test_without_a_normal_the_limitation_reads_word_for_word_as_before():
    out = global_limitations(_analysis(["tumour_vcf"]))
    assert LEVEL1_NORMAL_SENTENCE in out


def test_with_a_normal_the_limitation_says_what_the_subtraction_settles_and_what_it_does_not():
    out = global_limitations(_analysis(["tumour_vcf", "matched_normal"]))
    assert LEVEL1_NORMAL_SENTENCE not in out
    line = next(x for x in out if x.startswith("A matched normal genome was supplied"))
    assert "somatic by that subtraction" in line
    assert "The match is exact" in line


def _tv(gnomad_af: float | None = None):
    from genomeos.cancer.tumour import TumourVariant

    return TumourVariant("chr12", 25245350, "C", "A", gene="KRAS", gnomad_af=gnomad_af)


def test_an_origins_basis_names_the_normal_only_when_one_was_held():
    assert origin_of(_tv(), "s", {}, None, {}).somatic_basis == LEVEL1_BASIS
    held = origin_of(_tv(), "s", {}, None, {}, matched_normal=True).somatic_basis
    assert held.startswith("absent from the matched normal's calls")


# --- the normal outranks the population rule ---------------------------------------------


def _vcf(path: Path, rows: list[tuple[str, int, str, str]]) -> str:
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
    lines += [f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t0/1" for c, p, r, a in rows]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_grade_still_ends_its_germline_line_with_the_verdict_the_pipeline_replaces():
    from genomeos.genome.variants import Variant

    v = Variant("chr1", 999, "A", ("G",), ".", (0, 1), False, ".")
    vep = {"chr1:999:A:G": {"gene": "X", "consequence": "missense_variant", "gnomad_af": 0.2}}
    t = grade([v], vep, {})[0]
    assert t.likely_germline
    assert any(e.endswith(GERMLINE_VERDICT) for e in t.evidence)


def test_with_a_normal_a_common_variant_absent_from_it_is_somatic_and_one_present_is_removed(
    tmp_path, monkeypatch
):
    import genomeos.cancer.tumour as tumour

    common = ("chr1", 1000, "A", "G")  # gnomAD says inherited; the normal does not carry it
    shared = ("chr1", 2000, "C", "T")  # the normal carries it
    rare = ("chr1", 3000, "G", "A")
    t = _vcf(tmp_path / "t.vcf", [common, shared, rare])
    n = _vcf(tmp_path / "n.vcf", [shared])
    vep = {
        "chr1:999:A:G": {"gene": "AAA", "consequence": "missense_variant", "gnomad_af": 0.2},
        "chr1:1999:C:T": {"gene": "BBB", "consequence": "missense_variant", "gnomad_af": 0.3},
        "chr1:2999:G:A": {"gene": "CCC", "consequence": "synonymous_variant", "gnomad_af": None},
    }
    monkeypatch.setattr(tumour, "annotate_vep", lambda variants, log=None, append_to=None: vep)

    alone = pl.rank_variants(t, None, None, {})
    assert [v.gene for v in alone if v.likely_germline] == ["AAA", "BBB"]
    assert GERMLINE_AF <= 0.2

    paired = pl.rank_variants(t, n, None, {})
    assert {v.gene for v in paired} == {"AAA", "CCC"}  # BBB was in the normal
    aaa = next(v for v in paired if v.gene == "AAA")
    assert not aaa.likely_germline
    assert not any(e.endswith(GERMLINE_VERDICT) for e in aaa.evidence)
    assert aaa.evidence[0].startswith("population frequency 0.2 (gnomAD), but absent from the matched normal")
    # with nothing set aside the order is the score's, not "set-aside last"
    assert [v.score for v in paired] == sorted((v.score for v in paired), reverse=True)


# --- end to end, offline, from the demo's own VEP annotations ----------------------------

VEP_CACHE = DEMO / "vep_cache.jsonl"


def _cached_keys() -> set[str]:
    return {line.split("\t", 1)[0] for line in VEP_CACHE.read_text().splitlines() if line}


def test_the_demo_runs_offline_at_level_two():
    from genomeos.cancer.tumour import _key
    from genomeos.therapeutics import analyse_vcf

    # the demo carries its own VEP annotations, so this run never fetches: every tumour record is in it
    keys = _cached_keys()
    assert all(_key(v) in keys for v in iter_vcf(TUMOUR, None, pass_only=False))
    a = analyse_vcf(
        str(TUMOUR), normal_vcf=str(NORMAL), net=False, indirect=False, top_genes=4, vep_append_to=DEMO
    )
    assert a["data_level"]["level_reached"] == 2
    assert not any(m["input"] == "matched normal genome" for m in a["missing_data"])
    assert any("RNA-seq" in m["input"] for m in a["missing_data"])
    assert not any(c.tumour.expression.known for c in a["candidates"])
    assert all(
        o.somatic_basis.startswith("absent from the matched normal")
        for c in a["candidates"]
        for o in c.origins
    )
    limits = global_limitations(a)
    assert any(x.startswith("A matched normal genome was supplied") for x in limits)
    assert LEVEL1_NORMAL_SENTENCE not in limits
