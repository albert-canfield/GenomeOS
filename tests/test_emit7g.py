# SPDX-License-Identifier: AGPL-3.0-or-later
"""The census of what `compile.py` emits per element now and at the bar (genomeos/attribution/emit7g.py).

Every program and cache record here is synthetic, in the compiled generator's own line format and the
sweep's own record shape, so the counts are checked on the committed code they call
(readergate.rule_targets, targetchoice.chromosome_choice and candidate_list, targets.genes_at_bar)
without the generated store; the census script reads that store.
"""

from __future__ import annotations

import inspect

import pytest

from genomeos.attribution import compile as cp
from genomeos.attribution import emit7g as e7
from genomeos.attribution.direction_v2 import DIRECTION_RULE_DEFAULT, compiled_action
from genomeos.attribution.targets import genes_at_bar
from genomeos.predict.enhancer_target import MIN_EFFECT

PROGRAM = """module human.noncoding.chrT

gene AAA { symbol: AAA; evidence: curated "GENCODE v50" symbol only, a stub }
gene BBB { symbol: BBB; evidence: curated "GENCODE v50" symbol only, a stub }
region U_chrT_0 {
  locus: chrT:0-99
  role: structural, assembly gap, no sequence to attribute
  confidence: 1.00
}
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
  targets: AAA
  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect -0.3 log2 fold change on deletion, probability unavailable
}
rule E3 activates AAA { strength: 0.3; when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" effect -0.3 log2 fold change on deletion, probability unavailable }
element E1_measured {
  class: enhancer
  locus: chrT:100-200
  targets: AAA
  evidence: experimental "CRISPRi" links AAA in K562 from 1 pairs
  confidence: 0.90
}
rule E1_measured activates BBB { strength: 0.2; when: cell_type = K562; evidence: experimental "CRISPRi, silenced in K562, held-out split, evaluation only, never a feature"; confidence: 0.90 }
"""  # noqa: E501 - the generator's own line format, which runs past the limit

CODING = {"AAA", "BBB", "CCC", "DDD"}


def gene(name: str, drop: float, rise: float = 0.0) -> dict:
    return {
        "gene": name,
        "max_drop_log2fc": drop,
        "max_rise_log2fc": rise,
        "max_drop_tissue": "K562",
        "max_rise_tissue": "HepG2",
        "by_cell": {"K562": drop},
    }


#: E1: its compiled target AAA is the strongest; LNC1 is at the bar and not coding; CCC is at the bar on
#: a fall and DDD on a rise; FFF is under the bar. So E1 gains two targets. E2's compiled target BBB is
#: not the first coding gene at the bar (CCC is stronger), and it gains one. E3 is not cached.
RECORDS = {
    "E1": {
        "genes": [
            gene("AAA", -0.5),
            gene("LNC1", -0.45),
            gene("CCC", -0.3),
            gene("DDD", -0.02, 0.2),
            gene("FFF", -0.05),
        ]
    },
    "E2": {"genes": [gene("BBB", -0.4), gene("CCC", -0.6)]},
}


class FakeResponses:
    """`ElementResponses`' two readers on in-memory records, by the committed `genes_at_bar` rule."""

    def at_bar(self, chrom: str, element_id: str, min_effect: float = 0.1):
        return genes_at_bar(RECORDS.get(element_id), min_effect)

    def element(self, chrom: str, element_id: str):
        return RECORDS.get(element_id)


@pytest.fixture(scope="module")
def census() -> dict:
    return e7.chromosome_census("chrT", PROGRAM, CODING, FakeResponses())


def test_the_bar_is_the_compilers_own_object_and_is_never_copied() -> None:
    assert e7.BAR is MIN_EFFECT


def test_program_counts_splits_the_rule_lines_the_way_readergate_does() -> None:
    now = e7.program_counts(PROGRAM)
    assert now["rules_predicted"] == 3
    assert now["rules_experimental"] == 1
    assert now["rules_total"] == 4


def test_program_counts_states_every_block_kind_including_the_region_with_no_target() -> None:
    now = e7.program_counts(PROGRAM)
    assert now["element_blocks"] == 4
    assert now["region_blocks"] == 1
    assert now["gene_blocks"] == 2
    assert now["module_blocks"] == 1
    assert now["domain_blocks"] == 0


def test_the_experimental_element_block_is_counted_as_carrying_no_predicted_rule() -> None:
    now = e7.program_counts(PROGRAM)
    assert now["elements_with_a_predicted_rule"] == 3
    assert now["elements_without_a_predicted_rule"] == 1


def test_the_emission_adds_rules_and_the_count_is_the_candidate_lists_length(census: dict) -> None:
    # E1 -> AAA, CCC, DDD; E2 -> BBB, CCC; E3 -> AAA (not cached)
    assert census["rules_predicted_now"] == 3
    assert census["rules_predicted_after"] == 6
    assert census["rules_added"] == 3


def test_no_element_without_a_compiled_target_gains_one(census: dict) -> None:
    assert census["elements_gaining_a_target_that_had_none"] == 0
    assert census["elements_reached"] == 3


def test_the_elements_that_gain_and_keep_are_counted_and_add_up(census: dict) -> None:
    assert census["elements_gaining_a_target"] == 2
    assert census["elements_keeping_one_target"] == 1
    assert census["elements_gaining_a_target"] + census["elements_keeping_one_target"] == 3


def test_targets_per_element_after_is_a_histogram_of_the_candidate_lists(census: dict) -> None:
    assert census["targets_per_element_after"] == {"1": 1, "2": 1, "3": 1, "4-5": 0, "6-10": 0, "11+": 0}
    assert census["largest_candidate_list"] == 3


def test_every_choice_count_is_stated_even_at_zero(census: dict) -> None:
    for key in e7.CHOICE_KEYS:
        assert key in census["choice_counts"], key
    assert census["choice_counts"]["not_cached"] == 1
    assert census["choice_counts"]["compiled_target_not_first_at_the_bar"] == 1
    assert census["choice_counts"]["compiled_target_not_at_the_bar"] == 0


def test_a_candidate_row_carries_what_a_rule_line_needs() -> None:
    rows = e7.candidate_rules(FakeResponses(), "chrT", "E1", ["AAA"], CODING)
    assert [r["gene"] for r in rows] == ["AAA", "CCC", "DDD"]
    assert [r["head"] for r in rows] == [True, False, False]
    assert rows[1] == {
        "gene": "CCC",
        "rank": 1,
        "head": False,
        "log2": -0.3,
        "action": "activates",
        "tissue": "K562",
        "strength": 0.3,
    }
    # the rise side keeps its own track and reads as a repressed target
    assert rows[2]["action"] == "inhibits"
    assert rows[2]["tissue"] == "HepG2"
    assert rows[2]["log2"] == 0.2


def test_an_uncached_element_keeps_its_compiled_target_alone_with_no_cache_reading() -> None:
    rows = e7.candidate_rules(FakeResponses(), "chrT", "E3", ["AAA"], CODING)
    assert rows == [{"gene": "AAA", "rank": 0, "head": True}]


def test_combine_adds_every_count_and_recomputes_none(census: dict) -> None:
    total = e7.combine({"chrT": census, "chrU": census})
    assert total["chromosomes"] == 2
    assert total["rules_predicted_after"] == 2 * census["rules_predicted_after"]
    assert total["now"]["rules_total"] == 2 * census["now"]["rules_total"]
    assert total["targets_per_element_after"]["2"] == 2
    assert total["largest_candidate_list"] == census["largest_candidate_list"]


# ---- the two emissions `compile.py` can write (ROADMAP §5 item 7g) ----------------------------------


def test_the_default_emission_is_the_one_target_one_and_is_not_moved_here() -> None:
    assert cp.EMISSION_DEFAULT == cp.ONE_TARGET
    assert cp.EMISSIONS == (cp.ONE_TARGET, cp.EVERY_CODING_TARGET_AT_THE_BAR)
    # the default is also the signature's, so no caller gains the new emission by not passing one
    assert inspect.signature(cp.compile_chromosome).parameters["emission"].default == cp.ONE_TARGET


def test_element_targets_under_the_default_is_one_record_per_element_and_opens_no_file() -> None:
    elements = [
        {"id": "E1", "predicted_coding": {"gene": "AAA", "log2_fold_change": -0.5}},
        {"id": "E2", "predicted_coding": {"gene": "BBB", "log2_fold_change": -0.4}},
    ]
    got = cp.element_targets("chrT", elements, cp.ONE_TARGET)
    assert {k: [t["gene"] for t in v] for k, v in got.items()} == {"E1": ["AAA"], "E2": ["BBB"]}
    # identity, not a copy: the head rule line must be written from the element's own record
    assert got["E1"][0] is elements[0]["predicted_coding"]


def test_an_unknown_emission_is_refused_by_name() -> None:
    with pytest.raises(ValueError, match="no emission"):
        cp.element_targets("chrT", [], "every_target_anywhere")


def test_a_further_target_is_predicted_coding_shaped_so_every_helper_applies() -> None:
    """The point of the shape: the direction seam, the effect note and the cell label are the
    committed ones, so a further target's rule line differs from the compiled target's in numbers
    only."""
    row = {
        "gene": "CCC",
        "rank": 1,
        "head": False,
        "log2": -0.3,
        "action": "activates",
        "tissue": "K562",
        "strength": 0.3,
    }
    pc = cp.extra_target_pc(row)
    assert pc == {"gene": "CCC", "log2_fold_change": -0.3, "tissue": "K562", "action": "activates"}
    assert compiled_action(pc, DIRECTION_RULE_DEFAULT) == "activates"
    assert cp._effect_note(pc) == "effect -0.3 log2 fold change on deletion, probability unavailable"
    rise = cp.extra_target_pc({"gene": "DDD", "log2": 0.2, "action": "inhibits", "tissue": "HepG2"})
    assert compiled_action(rise, DIRECTION_RULE_DEFAULT) == "inhibits"
    assert rise["tissue"] == "HepG2"


def test_no_strength_word_is_invented_for_a_gene_beyond_the_head() -> None:
    """`basis` reads `pc.get('strength') or 'weak'`; the sweep states no strength word for a further
    target, so the record carries none rather than a word nothing measured."""
    pc = cp.extra_target_pc({"gene": "CCC", "log2": -0.3, "action": "activates", "tissue": "K562"})
    assert "strength" not in pc
