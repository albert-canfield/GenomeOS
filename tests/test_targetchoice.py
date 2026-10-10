# SPDX-License-Identifier: AGPL-3.0-or-later
"""The multi-target compiled choice (genomeos/attribution/targetchoice.py).

Every program and response cache here is synthetic, in the compiled generator's own line format and the
sweep's own record shape, so the choices are checked on the committed code they call
(signaltrace.predicted_text, links_of, ChromLookup; targets.genes_at_bar, ElementResponses.at_bar's
rule) without the generated store; the run reads that store.
"""

from __future__ import annotations

import math
import random

import pytest

from genomeos.attribution import compile as cp
from genomeos.attribution import crispri
from genomeos.attribution import signaltrace as st
from genomeos.attribution import targetchoice as tc
from genomeos.attribution.targets import genes_at_bar

PROGRAM = """module human.noncoding.chrT

gene AAA { symbol: AAA; evidence: curated "GENCODE v50" symbol only, a stub }
gene BBB { symbol: BBB; evidence: curated "GENCODE v50" symbol only, a stub }
gene CCC { symbol: CCC; evidence: curated "GENCODE v50" symbol only, a stub }
element E1 {
  class: enhancer
  locus: chrT:100-200
  targets: AAA
  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect -0.5 log2 fold change on deletion, probability unavailable
}
rule E1 activates AAA { strength: 0.5; when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" effect -0.5 log2 fold change on deletion, probability unavailable }
element E2 {
  class: enhancer
  locus: chrT:1000-1100
  targets: BBB
  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect -0.4 log2 fold change on deletion, probability unavailable
}
rule E2 activates BBB { strength: 0.4; when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" effect -0.4 log2 fold change on deletion, probability unavailable }
element E3 {
  class: enhancer
  locus: chrT:3000-3100
  targets: CCC
  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect -0.3 log2 fold change on deletion, probability unavailable
}
rule E3 activates CCC { strength: 0.3; when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" effect -0.3 log2 fold change on deletion, probability unavailable }
element E1_measured {
  class: enhancer
  locus: chrT:100-200
  targets: AAA
  evidence: experimental "CRISPRi" links AAA in K562 from 1 pairs
  confidence: 0.90
}
rule E1_measured activates BBB { strength: 0.2; when: cell_type = K562; evidence: experimental "CRISPRi, silenced in K562, held-out split, evaluation only, never a feature"; confidence: 0.90 }
"""  # noqa: E501 - the generator's own line format, which runs past the limit

CODING = {"AAA", "BBB", "CCC", "DDD", "EEE"}


def gene(name: str, drop: float, rise: float = 0.0) -> dict:
    return {"gene": name, "max_drop_log2fc": drop, "max_rise_log2fc": rise, "by_cell": {"K562": drop}}


#: the sweep's records. E1: its compiled target AAA is the strongest; LNC1 passes the bar and is not
#: coding; DDD passes on a fall and EEE on a rise; FFF is under the bar. E2: its compiled target BBB is
#: NOT the first coding gene at the bar (DDD is stronger). E3 is not cached.
RECORDS = {
    "E1": {
        "genes": [
            gene("AAA", -0.5),
            gene("LNC1", -0.45),
            gene("DDD", -0.3),
            gene("EEE", -0.02, 0.2),
            gene("FFF", -0.05),
        ]
    },
    "E2": {"genes": [gene("BBB", -0.4), gene("DDD", -0.6)]},
}


class FakeResponses:
    """ElementResponses' at_bar on in-memory records, by the same genes_at_bar rule."""

    def at_bar(self, chrom: str, element_id: str, min_effect: float = 0.1):
        return genes_at_bar(RECORDS.get(element_id), min_effect)


@pytest.fixture(scope="module")
def cl() -> tc.ChoiceLookup:
    return tc.chromosome_choice("chrT", PROGRAM, CODING, FakeResponses())


def pair(start: int, end: int, gene_: str, drop: float = 0.3, chrom: str = "chrT", regulated: bool = False):
    p = crispri.Pair(chrom, start, end, gene_, "K562", "test", 5000.0, 1.0, 1.0, regulated)
    p.features = {"deletion_drop": drop}
    return p


def test_candidates_are_the_compiled_target_then_coding_genes_at_the_bar_strongest_first(cl):
    assert cl.candidates["E1"] == ("AAA", "DDD", "EEE")  # LNC1 is not coding, FFF is under the bar
    assert cl.candidates["E2"] == ("BBB", "DDD")  # the compiled target first even where DDD is stronger
    assert cl.candidates["E3"] == ("CCC",)  # not cached: the compiled target alone


def test_the_bar_is_the_compilers_own_min_effect():
    from genomeos.predict.enhancer_target import MIN_EFFECT

    assert tc.MIN_EFFECT == MIN_EFFECT == 0.1
    assert tc.registered_fields()["min_effect"] == 0.1


def test_experimental_rules_are_removed_before_any_choice_reads_the_program(cl):
    assert cl.counts["experimental_rule_lines_removed"] == 1
    assert "E1_measured" not in cl.candidates  # its BBB rule is experimental and is never a candidate
    assert cl.counts["elements"] == 3 and cl.counts["compiled_targets_1"] == 3


def test_the_census_counts_name_each_silence(cl):
    assert cl.counts["not_cached"] == 1
    assert cl.counts["compiled_target_not_first_at_the_bar"] == 1  # E2: DDD is stronger than BBB
    assert cl.counts["compiled_target_not_at_the_bar"] == 0


def test_the_one_target_choice_is_signaltraces_target_choice(cl):
    look = st.ChromLookup("chrT")
    look.links = cl.look.links
    look.index({e: (s, t) for e, _c, s, t in cp.program_elements(PROGRAM) if e in cl.look.links})
    for start, end, g in [(150, 160, "AAA"), (150, 160, "DDD"), (1050, 1060, "BBB"), (5000, 5010, "AAA")]:
        p = pair(start, end, g)
        assert tc.trace_pair(p, cl)["one_target_kept"] == float(bool(look.named(start, end, g)))


def test_a_pair_on_another_compiled_gene_is_kept_only_by_the_wider_choices(cl):
    p = pair(150, 160, "DDD", drop=0.25)
    f = tc.trace_pair(p, cl, top_k=1)
    assert (f["one_target_kept"], f["multi_target_kept"], f["top_k_kept"]) == (0.0, 1.0, 0.0)
    assert f["multi_target_drop"] == 0.25 and f["one_target_drop"] == 0.0
    f2 = tc.trace_pair(p, cl, top_k=2)
    assert f2["top_k_kept"] == 1.0 and f2["top_k_drop"] == 0.25
    assert tc.rank_of(p, cl) == 2


def test_the_value_is_the_anchors_deletion_drop_on_kept_pairs_only(cl):
    for p in (pair(150, 160, "AAA", 0.4), pair(150, 160, "LNC1", 0.4), pair(150, 160, "FFF", 0.4)):
        f = tc.trace_pair(p, cl, top_k=3)
        for c in (tc.ONE, tc.MULTI, tc.TOPK):
            assert f[f"{c}_drop"] == (0.4 if f[f"{c}_kept"] else 0.0)
    assert tc.trace_pair(pair(150, 160, "LNC1"), cl)["multi_target_kept"] == 0.0  # not coding


def test_a_pair_on_no_compiled_element_is_kept_by_no_choice(cl):
    for c in (cl, None):
        f = tc.trace_pair(pair(9000, 9100, "AAA"), c, top_k=2)
        assert not any(f[f"{x}_kept"] for x in (tc.ONE, tc.MULTI, tc.TOPK))


def test_the_overlap_is_half_open(cl):
    assert tc.trace_pair(pair(200, 300, "AAA"), cl)["one_target_kept"] == 0.0
    assert tc.trace_pair(pair(199, 300, "AAA"), cl)["one_target_kept"] == 1.0


def test_the_harness_refuses_a_choice_that_drops_a_pair_the_narrower_one_keeps():
    p = pair(1, 2, "X")
    p.features.update({"one_target_kept": 1.0, "top_k_kept": 0.0, "multi_target_kept": 1.0})
    assert tc.harness([p], (tc.ONE, tc.MULTI, tc.TOPK)) == [
        "1 pairs kept by one_target are not kept by top_k"
    ]
    p.features["top_k_kept"] = 1.0
    assert tc.harness([p], (tc.ONE, tc.MULTI, tc.TOPK)) == []
    p.features["multi_target_kept"] = 0.0
    assert tc.harness([p], (tc.ONE, tc.MULTI)) == ["1 pairs kept by one_target are not kept by multi_target"]


@pytest.mark.parametrize(
    ("added", "vs_one", "vs_anchor", "outcome"),
    [
        (0, [0.01, 0.2], [-0.1, -0.01], "NONE_ADDED"),
        (5, [0.01, 0.2], [-0.1, -0.01], "BETTER_BELOW_ANCHOR"),
        (5, [0.01, 0.2], [-0.1, 0.05], "BETTER_NO_DIFFERENCE_FROM_ANCHOR"),
        (5, [0.01, 0.2], [0.01, 0.05], "BETTER_ABOVE_ANCHOR"),
        (5, [-0.01, 0.2], [-0.1, -0.01], "NO_BETTER"),
        (5, [0.0, 0.2], [-0.1, -0.01], "NO_BETTER"),
        (5, [-0.2, -0.01], [-0.3, -0.1], "WORSE"),
    ],
)
def test_the_decision_follows_the_intervals(added, vs_one, vs_anchor, outcome):
    assert tc.decide(added, vs_one, vs_anchor) == outcome
    assert (outcome in tc.FALSIFIED_BY) == (outcome not in tc.BETTER)


def test_every_outcome_sentence_is_filled_by_the_numbers_and_nothing_else():
    interval = {"gain": 0.05, "ci95": [0.01, 0.1]}
    row = {
        "population": {"pairs": 1918, "kept": 400, "kept_regulated": 80},
        "gain_over_activity_distance": interval,
        "against_one_target": interval,
        "against_anchor": {"gain": -0.02, "ci95": [-0.05, 0.01]},
    }
    one = {"population": {"kept": 39, "kept_regulated": 35}, "gain_over_activity_distance": interval}
    nums = tc.numbers(row, one, tc.phrase_of(tc.MULTI, None))
    for kind in tc.OUTCOMES:
        assert "{" not in tc.reading(kind, what="x", **nums)
    text = tc.reading("NO_BETTER", **nums)
    assert "1,918" in text and "+0.0500 [95% CI +0.0100, +0.1000]" in text and "-0.0200" in text
    assert tc.phrase_of(tc.TOPK, 3) == "keeping the first 3 candidate targets per element"


def test_the_registered_fields_round_trip_and_bind_every_constant():
    f = tc.registered_fields()
    assert set(f) == set(tc.BOUND)
    assert f["anchor"] == st.ANCHOR and f["bootstrap"] == st.BOOTSTRAP
    assert f["base_columns"] == list(crispri.FEATURES["activity + distance"])
    assert f["comparator"]["gain"] == 0.0792 and f["comparator"]["kept"] == 39


def test_k_is_the_smallest_with_the_largest_in_sample_average_precision():
    by_k = {str(k): {"average_precision_in_sample": 0.5} for k in tc.K_TRIED}
    by_k["3"]["average_precision_in_sample"] = by_k["5"]["average_precision_in_sample"] = 0.6
    assert tc.choose_k(by_k) == 3
    assert tc.secondary_status(1, 9).startswith("not read: TOP_K is 1")
    assert tc.secondary_status(3, 3).startswith("not read: no compiled element")
    assert tc.secondary_status(3, 4) == "read beside the primary"
    assert tc.read_choices(None) == tc.read_choices(1) == (tc.ONE, tc.MULTI)
    assert tc.read_choices(3) == (tc.ONE, tc.MULTI, tc.TOPK)


def test_the_comparator_check_names_every_pair_that_differs_from_m1c(cl):
    p, q = pair(150, 160, "AAA", 0.4), pair(150, 160, "DDD", 0.4)
    for x in (p, q):
        x.features.update(tc.trace_pair(x, cl))
    table = {
        "columns": ["pair", "regulated", "s2_target_kept", "s2_target_drop"],
        "rows": [[st.pair_id(p), False, 1.0, 0.4], [st.pair_id(q), False, 1.0, 0.4]],
    }
    assert tc.comparator_mismatches([p, q], table) == [st.pair_id(q)]


def population(n_chroms: int = 12, per: int = 40, seed: int = 1) -> list[crispri.Pair]:
    """Synthetic pairs over enough chromosomes for gain_interval, with the base columns filled."""
    rng = random.Random(seed)
    out = []
    for c in range(n_chroms):
        for i in range(per):
            regulated = rng.random() < 0.15
            p = pair(i * 10, i * 10 + 5, f"G{i}", chrom=f"chr{c + 1}", regulated=regulated)
            d = rng.uniform(1_000, 500_000)
            p.features = {
                "log_distance": math.log(d),
                "log_activity": rng.random(),
                "activity_over_distance": rng.random() - math.log(d),
                "top_target": 0.0,
                "deletion_drop": rng.random() * (2 if regulated else 1),
            }
            one = regulated and rng.random() < 0.2
            multi = one or (rng.random() < (0.6 if regulated else 0.2))
            for choice, k in ((tc.ONE, one), (tc.MULTI, multi)):
                p.features[f"{choice}_kept"] = float(k)
                p.features[f"{choice}_drop"] = p.features["deletion_drop"] if k else 0.0
            out.append(p)
    return out


def test_the_score_reads_every_choice_against_the_one_target_and_the_anchor():
    fit, pop = population(seed=1), population(seed=2)
    _, base = st.fit_scores(fit, pop, tc.BASE_COLUMNS)
    w, anchor_scores = st.fit_scores(fit, pop, tc.ANCHOR_COLUMNS)
    anchor = {"base": base, "scores": anchor_scores, "weights": w}
    one = tc.score_one(fit, pop, anchor)
    out = tc.score(fit, pop, anchor, one, (tc.ONE, tc.MULTI), None)
    row = out["choices"][tc.MULTI]
    assert set(out["choices"]) == {tc.MULTI}
    assert row["outcome"] == out["outcome"] in tc.OUTCOMES
    assert row["population_pairs_added_over_one_target"] > 0
    assert out["falsified"] == (out["outcome"] in tc.FALSIFIED_BY)
    assert "{" not in out["reading"]
    assert not tc.comparator_reproduced(one)  # synthetic numbers are not M1c's


def test_a_multi_choice_that_adds_nothing_is_none_added():
    fit, pop = population(seed=3), population(seed=4)
    for p in fit + pop:
        p.features["multi_target_kept"] = p.features["one_target_kept"]
        p.features["multi_target_drop"] = p.features["one_target_drop"]
    _, base = st.fit_scores(fit, pop, tc.BASE_COLUMNS)
    anchor = {"base": base, "scores": base}
    one = tc.score_one(fit, pop, anchor)
    out = tc.score(fit, pop, anchor, one, (tc.ONE, tc.MULTI), None)
    assert out["outcome"] == "NONE_ADDED" and out["falsified"]


def test_top_k_is_the_k_the_committed_training_census_names():
    import json
    from pathlib import Path

    import tracked_paths as tp

    tp.must_be_present(tc.CENSUS)
    census = json.loads(Path(tc.CENSUS).read_text())
    assert census["top_k"] == tc.TOP_K == 6
    assert census["heldout_read"] is False
    assert census["census_rule"] == tc.CENSUS_RULE
    assert tc.read_choices(tc.TOP_K) == (tc.ONE, tc.MULTI, tc.TOPK)


def test_the_run_pins_its_own_code_and_reads_the_population_by_m1cs_ids():
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import targetchoice_run as run

    assert set(run.FROZEN) == {
        "genomeos/attribution/targetchoice.py",
        "scripts/targetchoice_census.py",
        "scripts/targetchoice_register.py",
        "scripts/targetchoice_run.py",
    }
    assert "genomeos/attribution/signaltrace.py" in run.STAGE_CODE
    training = [pair(1, 2, "A", chrom="chr2"), pair(1, 2, "B", chrom="chr10")]
    assert run.chromosomes_of(training, ["chrX:5-9:G", "chr2:1-2:A"]) == ["chr2", "chr10", "chrX"]
