# SPDX-License-Identifier: AGPL-3.0-or-later
"""The second-targets test against Gasperini 2019 (genomeos/attribution/secondtargets.py), on synthetic
tables and a synthetic compiled program: the pair table's filters, the independence exclusion, the arms,
the matching, the readout's hit rule, the decision and every outcome sentence."""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import pytest

from genomeos.attribution import secondtargets as s2
from genomeos.attribution import targetchoice as tc
from genomeos.attribution.targets import genes_at_bar

ROOT = Path(__file__).resolve().parents[1]

PROGRAM = """module human.noncoding.chrT

gene AAA { symbol: AAA; evidence: curated "GENCODE v50" symbol only, a stub }
gene BBB { symbol: BBB; evidence: curated "GENCODE v50" symbol only, a stub }
element E1 {
  class: enhancer
  locus: chrT:1000-1200
  targets: AAA
  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect -0.5 log2 fold change on deletion, probability unavailable
}
rule E1 activates AAA { strength: 0.5; when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" effect -0.5 log2 fold change on deletion, probability unavailable }
element E2 {
  class: enhancer
  locus: chrT:5000-5100
  targets: BBB
  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect -0.4 log2 fold change on deletion, probability unavailable
}
rule E2 activates BBB { strength: 0.4; when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" effect -0.4 log2 fold change on deletion, probability unavailable }
"""  # noqa: E501 - the generator's own line format

CODING = {"AAA", "BBB", "CCC", "DDD", "EEE", "FFF", "GGG"}


def gene(name: str, drop: float) -> dict:
    return {"gene": name, "max_drop_log2fc": drop, "max_rise_log2fc": 0.0, "by_cell": {"K562": drop}}


#: E1's compiled target AAA, then CCC and DDD at the bar (CCC stronger), EEE under it
RECORDS = {"E1": {"genes": [gene("AAA", -0.5), gene("DDD", -0.2), gene("CCC", -0.3), gene("EEE", -0.05)]}}


class FakeResponses:
    def at_bar(self, chrom: str, element_id: str, min_effect: float = 0.1):
        return genes_at_bar(RECORDS.get(element_id), min_effect)


@pytest.fixture(scope="module")
def cl() -> tc.ChoiceLookup:
    return tc.chromosome_choice("chrT", PROGRAM, CODING, FakeResponses())


HEADER = (
    "gRNAgroup.chr\tgRNAgroup.start\tgRNAgroup.stop\tgRNAgroup\tgeneral_group\tchr.targetgene\t"
    "start.targetgene\tstop.targetgene\tENSG.targetgene\ttargetgene_short_name\tstrand.targetgene\tpairs"
)


def line(chrom, start, stop, group, sel, tss, ensg, sym) -> str:
    return "\t".join(
        [chrom, str(start), str(stop), group, sel, "chrT", str(tss), str(tss + 1), ensg, sym, "+", "x"]
    )


def identity(chrom: str, start: int, end: int):
    return None if chrom == "chrLOST" else (chrom, start, end)


def screen_test(group, ensg, symbol, distance, start=1000, end=1200, chrom="chrT") -> s2.Tested:
    return s2.Tested(group, chrom, start, end, ensg, symbol, distance)


# ---- the pair table


def test_the_pair_table_keeps_top_two_candidate_tests_and_counts_every_exclusion():
    rows = [
        HEADER,
        line("NTC", "NTC", "NTC", "bassik_mch", "NTC", 10, "ENSG1", "A1"),
        line("chrT", 1000, 1200, "chrT.1_top_two", "test_5k-untested", 21100, "ENSG00000000001.5", "AAA"),
        line("chrT", 1000, 1200, "chrT.1_top_two", "test_5k-untested", 21100, "ENSG00000000001.5", "AAA"),
        line("chrT", 1000, 1200, "chrT.1_second_two", "test_5k-untested", 21100, "ENSG00000000001", "AAA"),
        line("chrT", 900, 901, "AAA_TSS", "TSS", 900, "ENSG00000000001", "AAA"),
        line("chrT", 900, 901, "AAA_TSS", "TSS", 5000, "ENSG00000000002", "BBB"),
        line("chrT", 1, 2, "chrT.9_control_HS2", "positive_ctrl", 50, "ENSG9", "HBG1"),
        line("chrT", 3000, 3200, "chrT.2_top_two", "not_model_selected_K1119hit-hit", 3100, "ENSG3", "CCC"),
        line("chrLOST", 3000, 3200, "chrLOST.3_top_two", "test_5k-hit", 3100, "ENSG4", "DDD"),
    ]
    tests, tss, counts = s2.read_pair_table([r + "\n" for r in rows], identity)
    assert [t.key for t in tests] == ["chrT.1_top_two:ENSG00000000001"]
    assert tests[0].distance == 20000  # |21100 - 1100|, in the table's own build
    assert [t.key for t in tss] == ["AAA_TSS:ENSG00000000001"]  # the TSS group's own gene only
    assert counts["non_targeting"] == 1 and counts["other_guide_group"] == 1 and counts["other_control"] == 1
    assert counts["tss_control"] == 2 and counts["selected_for_a_pilot_hit"] == 1
    assert counts["not_lifted"] == 1 and counts["repeated"] == 1 and counts["tests"] == 1
    assert counts["candidate_enhancer_tests"] == 4


def test_the_pair_table_refuses_a_table_without_its_columns():
    with pytest.raises(KeyError):
        s2.read_pair_table(["gRNAgroup.chr\tgRNAgroup.start\n"], identity)


def test_ensembl_ids_lose_prefix_and_version():
    assert s2.ensg_of("ENSEMBL:ENSG00000108179.12") == "ENSG00000108179"
    assert s2.ensg_of("ENSG00000108179") == "ENSG00000108179"
    assert s2.ensg_of("NA") == ""


# ---- independence

BENCH = [
    "chrom\tchromStart\tchromEnd\tmeasuredGeneSymbol\tmeasuredGeneEnsemblId\tDataset\tRegulated\n",
    "chrT\t1100\t1300\tAAA\tENSG00000000001\tGasperini2019\tTRUE\n",
    "chrT\t1150\t1160\tOLDNAME\tENSG00000000005\tNasser2021\tFALSE\n",
    "chrT\t1190\t1250\tZZZ\tENSG00000000099\tXie\tFALSE\n",
]


def test_a_test_is_excluded_when_the_benchmark_names_its_gene_at_an_overlapping_element():
    index = s2.BenchmarkIndex.from_tables([BENCH])
    assert s2.independence(screen_test("g", "ENSG00000000001", "AAA", 10), index) == (
        "same_gene",
        ["Gasperini2019"],
    )
    # the same gene by Ensembl id under another symbol
    assert s2.independence(screen_test("g", "ENSG00000000005", "NEW", 10), index)[0] == "same_gene"
    # the same element, another gene: independent pair by pair, not element by element
    assert s2.independence(screen_test("g", "ENSG00000000002", "BBB", 10), index) == (
        "element_other_gene",
        [],
    )
    # half-open: [1300, 1400) does not touch [1100, 1300)
    far = screen_test("g", "ENSG00000000001", "AAA", 10, start=1300, end=1400)
    assert s2.independence(far, index) == ("element_untested", [])
    assert index.datasets == Counter({"Gasperini2019": 1, "Nasser2021": 1, "Xie": 1})


def test_the_benchmark_is_read_without_its_label_columns():
    assert not {"Regulated", "Significant", "EffectSize", "pValueAdjusted"} & set(s2.BENCHMARK_COLUMNS)
    assert "Dataset" in s2.BENCHMARK_COLUMNS


# ---- arms and matching


def test_the_arms_name_the_compiled_target_the_shortlist_and_the_rest(cl):
    assert cl.candidates["E1"] == ("AAA", "CCC", "DDD")  # compiled target, then the bar, strongest first
    t = screen_test("g", "E", "X", 10)
    assert s2.arm_of(t, "AAA", cl) == (s2.COMPILED, 1)
    assert s2.arm_of(t, "CCC", cl) == (s2.SHORT, 2)
    assert s2.arm_of(t, "DDD", cl) == (s2.SHORT, 3)
    assert s2.arm_of(t, "EEE", cl) == (s2.OTHER, None)  # under the bar
    off = screen_test("g", "E", "X", 10, start=2000, end=2100)
    assert s2.arm_of(off, "CCC", cl) == (None, None)
    assert s2.arm_of(t, "CCC", None) == (None, None)


def test_the_shortlist_census_counts_second_and_further_targets(cl):
    out = s2.shortlist_census({"chrT": cl})
    assert out["compiled_elements"] == 2
    assert out["with_a_second_target"] == 1 and out["second_and_further_targets"] == 2


def test_rows_drop_benchmark_tests_non_coding_genes_and_elements_without_a_compiled_element(cl):
    index = s2.BenchmarkIndex.from_tables([BENCH])
    coding = {"chrT": {"ENSG00000000001": "AAA", "ENSG00000000003": "CCC", "ENSG00000000007": "GGG"}}
    tests = [
        screen_test("g1", "ENSG00000000001", "AAA", 100),  # in the benchmark
        screen_test("g1", "ENSG00000000003", "CCC", 100),  # shortlist
        screen_test("g1", "ENSG00000000007", "GGG", 150),  # other
        screen_test("g1", "ENSG00000000042", "LNC", 150),  # not coding
        screen_test("g2", "ENSG00000000003", "CCC", 100, start=7000, end=7100),  # no compiled element
    ]
    rows, split = s2.rows_of(tests, index, coding, {"chrT": cl})
    assert [(r.symbol, r.arm) for r in rows] == [("CCC", s2.SHORT), ("GGG", s2.OTHER)]
    assert split["independence"] == {"same_gene": 1, "element_other_gene": 3, "element_untested": 1}
    assert split["after_independence"]["not_protein_coding_in_gencode_v50"] == 1
    assert split["after_independence"]["no_compiled_element_overlaps"] == 1
    assert split["excluded_by_benchmark_dataset"] == {"Gasperini2019": 1}
    assert not any(r.element_untested for r in rows)


def row(group, ensg, arm, distance, rank=None, chrom="chrT") -> s2.Row:
    return s2.Row(screen_test(group, ensg, ensg, distance, chrom=chrom), ensg, arm, rank, False)


def test_matching_takes_every_other_gene_of_the_element_within_a_factor_of_two():
    rows = [
        row("g1", "S1", s2.SHORT, 10_000, 2),
        row("g1", "O1", s2.OTHER, 19_000),  # within x2
        row("g1", "O2", s2.OTHER, 5_000),  # exactly half: within
        row("g1", "O3", s2.OTHER, 21_000),  # beyond x2
        row("g1", "C1", s2.COMPILED, 3_000, 1),  # nearest gene of the element
        row("g2", "S2", s2.SHORT, 400_000, 3),
        row("g2", "O4", s2.OTHER, 100_000),  # beyond x2: S2 is unmatched
        row("g3", "O5", s2.OTHER, 10_000),  # an element without a shortlist test adds nothing
    ]
    units, unmatched = s2.matched(rows)
    assert unmatched == 1 and len(units) == 1
    u = units[0]
    assert u.test.symbol == "S1"
    assert {c.symbol for c in u.controls} == {"O1", "O2"}
    assert u.nearest.symbol == "C1" and [c.symbol for c in u.compiled] == ["C1"]
    assert s2.population_ids(units) == ["g1:S1"]
    assert len(s2.ids_sha256(["g1:S1"])) == 64


def test_the_census_counts_the_population_label_free():
    rows = [row("g1", "S1", s2.SHORT, 10_000, 2), row("g1", "O1", s2.OTHER, 12_000)]
    units, unmatched = s2.matched(rows)
    split = {"independence": {"same_gene": 3, "element_other_gene": 1, "element_untested": 1}}
    c = s2.census(Counter({"tests": 5}), rows, split, units, unmatched)
    assert c["excluded_same_gene_in_benchmark"] == 3 and c["independent"] == 2
    assert c["shortlist_tests"] == 1 and c["matched_shortlist_tests"] == 1
    assert c["shortlist_by_rank"] == {"2": 1}
    assert json.loads(json.dumps(c)) == c  # the registration stores it as JSON and the run compares


# ---- the readout

READOUT = [
    "pairs4merge\tbeta\tfold_change.transcript_remap\tpvalue.empirical\tpvalue.empirical.adjusted\n",
    "g1:S1\t-1\t0.6\t0.001\t0.01\n",  # hit
    "g1:O1\t1\t1.4\t0.001\t0.01\n",  # significant increase: not a hit
    "g1:O2\t0\t0.9\t0.5\t0.6\n",  # not a hit
    "g1:C1\t0\tNA\tNA\tNA\n",  # not measured
    "g9:X\t-1\t0.5\t0.001\t0.01\n",  # not wanted, not read into the result
]


def test_the_hit_rule_is_the_screens_own_call_a_decrease():
    hits = s2.read_readout(READOUT, {"g1:S1", "g1:O1", "g1:O2", "g1:C1"})
    assert hits == {"g1:S1": True, "g1:O1": False, "g1:O2": False, "g1:C1": None}


def test_a_readout_without_the_registered_columns_stops_before_any_value():
    with pytest.raises(KeyError):
        s2.read_readout(["pairs4merge\tpvalue\tfold_change\n", "g1:S1\t0.01\t0.5\n"], {"g1:S1"})


def test_the_positive_control_is_the_share_of_measured_tss_groups_that_are_hits():
    tss = [
        screen_test("A_TSS", "EA", "A", 1),
        screen_test("B_TSS", "EB", "B", 1),
        screen_test("C_TSS", "EC", "C", 1),
    ]
    out = s2.positive_control(tss, {"A_TSS:EA": True, "B_TSS:EB": False, "C_TSS:EC": None})
    assert out == {"groups": 3, "measured": 2, "hits": 1, "share": 0.5}


def test_unmeasured_tests_leave_the_comparison_and_are_counted():
    rows = [
        row("g1", "S1", s2.SHORT, 10_000, 2),
        row("g1", "O1", s2.OTHER, 12_000),
        row("g1", "O2", s2.OTHER, 9_000),
        row("g2", "S2", s2.SHORT, 10_000, 2),
        row("g2", "O3", s2.OTHER, 10_000),
    ]
    units, _ = s2.matched(rows)
    hits = {"g1:S1": True, "g1:O1": False, "g1:O2": None, "g2:S2": None, "g2:O3": True}
    scored, lost = s2.scored_units(units, hits)
    assert len(scored) == 1 and scored[0].hit and scored[0].control_rate == 0.0
    assert lost == Counter({"shortlist_test_not_measured": 1, "controls_not_measured": 1})


def test_the_interval_resamples_whole_chromosomes_and_is_reproducible():
    units = [s2.Scored(f"chr{i % 5}", i % 3 == 0, 0.1, False, False, 2) for i in range(60)]
    a = s2.chromosome_interval(units, lambda u: u.control_rate, n=200)
    b = s2.chromosome_interval(units, lambda u: u.control_rate, n=200)
    assert a == b and a["chromosomes"] == 5
    assert a["ci95"][0] <= a["difference"] <= a["ci95"][1]
    assert a["difference"] == pytest.approx(100 * (20 - 6) / 60, abs=1e-3)
    assert s2.chromosome_interval([], lambda u: 0.0)["ci95"] is None


@pytest.mark.parametrize(
    "matched_n, held, ci, outcome",
    [
        (99, True, (1.0, 2.0), "TOO_FEW_PAIRS"),
        (99, False, None, "TOO_FEW_PAIRS"),
        (100, False, None, "NOT_RUN_NO_READOUT"),
        (100, True, (0.1, 2.0), "MORE_OFTEN"),
        (100, True, (-0.1, 2.0), "NO_DIFFERENCE"),
        (100, True, (0.0, 2.0), "NO_DIFFERENCE"),
        (100, True, (-2.0, -0.1), "LESS_OFTEN"),
        (100, True, None, "HARNESS_DEFECT"),
    ],
)
def test_the_decision_follows_the_registered_order(matched_n, held, ci, outcome):
    assert s2.decide(matched_n, held, ci) == outcome
    assert set(s2.FALSIFIED_BY) <= set(s2.RAN) and "MORE_OFTEN" not in s2.FALSIFIED_BY


def test_every_outcome_sentence_is_filled_by_its_numbers_and_nothing_else():
    census = {
        "candidate_enhancer_tests_read": 128737,
        "excluded_same_gene_in_benchmark": 6015,
        "independent": 122722,
        "coding_tests_at_compiled_elements": 57112,
        "shortlist_tests": 1411,
        "matched_shortlist_tests": 872,
    }
    nums = s2.census_numbers(census)
    for kind in ("TOO_FEW_PAIRS", "NOT_RUN_NO_READOUT"):
        text = s2.reading(kind, **nums)
        assert "{" not in text and "872" in text and "6,015" in text
    assert "0 of the 872" in s2.reading("NOT_RUN_NO_READOUT", **nums)
    ran = dict(n="872", h=30, r="3.44", rc="1.10", d="+2.34", lo="+0.90", hi="+3.80", rn="9.00", dn="-5.56")
    ran.update(nlo="-8.00", nhi="-3.00")
    for kind in s2.RAN:
        assert "{" not in s2.reading(kind, **ran)
    assert "not function elsewhere" in s2.reading("MORE_OFTEN", **ran)
    assert "is not no difference" in s2.reading("NO_DIFFERENCE", **ran)
    assert "{" not in s2.reading("HARNESS_DEFECT", what="x")


def test_the_registered_fields_round_trip_and_bind_every_constant():
    fields = s2.registered_fields()
    assert set(fields) == set(s2.BOUND)
    assert json.loads(json.dumps(fields)) == fields
    assert fields["min_effect"] == tc.MIN_EFFECT == 0.1
    assert fields["hit"]["p_below"] == 0.1 and fields["distance_factor"] == 2.0


def test_the_run_pins_its_own_code_and_the_code_it_calls():
    sys.path.insert(0, str(ROOT / "scripts"))
    import secondtargets_run as run

    assert set(run.FROZEN) == {
        "genomeos/attribution/secondtargets.py",
        "scripts/secondtargets_register.py",
        "scripts/secondtargets_run.py",
    }
    assert {"genomeos/attribution/targetchoice.py", "genomeos/attribution/measured_unknown.py"} <= set(
        run.STAGE_CODE
    )


def test_the_committed_registration_binds_this_module():
    import tracked_paths as tp

    path = "data/results/secondtargets_registration.json"
    tp.must_be_present(path)
    reg = json.loads((ROOT / path).read_text())
    assert {f: reg[f] for f in s2.BOUND} == s2.registered_fields()
    assert reg["readout_values_read"] == 0 and reg["readout"]["held"] is False
    assert s2.ids_sha256(reg["population_ids"]) == reg["population_sha256"]
    assert reg["census"]["matched_shortlist_tests"] == len(reg["population_ids"]) == 872
    assert reg["census"]["excluded_same_gene_in_benchmark"] + reg["census"]["independent"] == 128737
    assert reg["frozen_code"]["commit"].startswith("d365cdd")
