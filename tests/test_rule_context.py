# SPDX-License-Identifier: AGPL-3.0-or-later
"""R1 of the external review of 2026-09-28: context in executable rules.

The compiler wrote the cell a screen used only into evidence text, so a K562 rule ran in every
simulated cell, and the measured layer kept one rule per element and gene, the strongest across
cells. These tests pin the fix pre-registered in docs/ATTRIBUTION.md (R1, 2026-09-28): one rule per
(element, gene, cell), each gated by an executable `when: cell_type = <cell>`, and a rule with no
recorded context gated on `unknown`, which matches no cell.
"""

from __future__ import annotations

import json
from pathlib import Path

from genomeos.attribution import compile as cp
from genomeos.attribution import measured
from genomeos.ir import Module
from genomeos.lang import parse
from genomeos.runtime.grn import NetworkRuntime

CHROM = "chrT"


def pair(gene, effect, cell, split=measured.TRAINING, significant=True, power=0.95):
    return measured.CrispriPair(
        chrom=CHROM,
        start=1000,
        end=1300,
        gene=gene,
        cell=cell,
        dataset="Syn2026",
        reference="synthetic",
        regulated=significant and effect < 0,
        significant=significant,
        effect_size=effect,
        p_adjusted=0.01 if significant else 0.5,
        split=split,
        power_at_effect_size_20=power,
    )


def row_of(*pairs):
    m = measured.Layer(chrom=CHROM, crispri=list(pairs)).for_element(1000, 1300)
    return {
        "id": "E1",
        "start": 1000,
        "end": 1300,
        "class": "enhancer",
        "measured": m,
        "agreement": measured.agreement("AAA", m),
        "predicted_gene": "AAA",
        "confidence": measured.confidence_of(m),
    }


def program(row, genes=("AAA",)):
    ident_of = {g: cp.ident(g) for g in genes}
    block = cp._measured_blocks(CHROM, row, {}, ident_of)
    stubs = [
        f'gene {i} {{ symbol: {g}; evidence: curated "x"; confidence: 0.9 }}' for g, i in ident_of.items()
    ]
    return "\n".join(["module t", *stubs, *block]) + "\n", block


def measured_rules(module):
    return [r for r in module.rules if r.source == "E1_measured"]


# --- the review's acceptance tests, as it words them ------------------------------------------------
def test_a_k562_specific_rule_is_inactive_in_hepg2():
    text, _ = program(row_of(pair("AAA", -0.3, "K562")))
    module = parse(text)
    (rule,) = measured_rules(module)
    assert rule.when == {measured.CONTEXT_KEY: "K562"}
    assert not [r for r in module.active_rules({"cell_type": "HepG2"}) if r.source == "E1_measured"]
    assert [r for r in module.active_rules({"cell_type": "K562"}) if r.source == "E1_measured"] == [rule]
    # and in the runtime that integrates rules, not only in the filter it calls
    assert not NetworkRuntime(module, context={"cell_type": "HepG2"}).active_rules
    assert NetworkRuntime(module, context={"cell_type": "K562"}).active_rules == [rule]


def round_trip(text: str) -> Module:
    """Compiled text, parse, BioIR dict, JSON, and back."""
    return Module.from_dict(json.loads(json.dumps(parse(text).to_dict())))


def test_conflicting_results_from_two_cell_types_survive_compilation_and_round_trip():
    text, _ = program(row_of(pair("AAA", -0.2, "K562"), pair("AAA", -0.7, "HepG2")))
    for module in (parse(text), round_trip(text)):
        rules = {r.when["cell_type"]: r for r in measured_rules(module)}
        assert set(rules) == {"K562", "HepG2"}  # two cells, two rules, never the strongest of them
        assert (rules["K562"].strength, rules["HepG2"].strength) == (0.2, 0.7)
        assert not [r for r in module.active_rules({"cell_type": "HepG2"}) if r is rules["K562"]]


def test_a_null_in_one_cell_beside_a_link_in_another_is_kept_cell_by_cell():
    """A well-powered null in HepG2 makes no rule, but it is no longer erased by K562's link."""
    text, block = program(row_of(pair("AAA", -0.3, "K562"), pair("AAA", 0.0, "HepG2", significant=False)))
    module = round_trip(text)
    assert [r.when for r in measured_rules(module)] == [{"cell_type": "K562"}]
    assert not [r for r in module.active_rules({"cell_type": "HepG2"}) if r.source == "E1_measured"]
    basis = next(line for line in block if "basis:" in line)
    assert "in HepG2 measured no effect on AAA" in basis


def test_a_significant_increase_in_another_cell_is_named_and_still_raises_no_rule():
    text, block = program(row_of(pair("AAA", -0.3, "K562"), pair("AAA", 0.4, "HepG2")))
    assert [r.when for r in measured_rules(parse(text))] == [{"cell_type": "K562"}]
    assert "in HepG2 significantly increased AAA" in next(line for line in block if "basis:" in line)


# --- a missing context is explicit ------------------------------------------------------------------
def test_a_rule_with_no_recorded_cell_is_gated_on_unknown_and_matches_no_cell():
    text, _ = program(row_of(pair("AAA", -0.3, "")))
    (rule,) = measured_rules(parse(text))
    assert rule.when == {"cell_type": measured.CONTEXT_UNKNOWN}
    for context in ({}, {"cell_type": "K562"}, {"cell_type": "unknown"}, {"cell_type": "any"}):
        assert not rule.applies(context), context


def test_every_compiled_rule_names_its_context(tmp_path):
    """Predicted rules too: the AlphaGenome track's biosample, or `unknown` where none was recorded."""
    elements = [
        {
            "id": "E1",
            "start": 1000,
            "end": 1300,
            "predicted_coding": {
                "gene": "AAA",
                "log2_fold_change": -0.5,
                "action": "activates",
                "tissue": "K562",
            },
        },
        {
            "id": "E2",
            "start": 5000,
            "end": 5300,
            "predicted_coding": {
                "gene": "BBB",
                "log2_fold_change": -0.4,
                "action": "activates",
                "tissue": "CD8-positive, alpha-beta T cell",
            },
        },
        {
            "id": "E3",
            "start": 9000,
            "end": 9300,
            "predicted_coding": {"gene": "CCC", "log2_fold_change": 0.3},
        },
    ]
    d = tmp_path / "results"
    d.mkdir()
    (d / f"budget_{CHROM}.json").write_text(json.dumps({"unknown_bp": 0, "blocks": []}))
    (d / f"enhancer_targets_{CHROM}.json").write_text(json.dumps({"elements": elements}))
    layer = measured.Layer(chrom=CHROM, crispri=[pair("AAA", -0.3, "K562"), pair("AAA", -0.6, "WTC11")])
    text = cp.compile_chromosome(CHROM, d, layer=layer)
    module = parse(text)
    assert all(r.when.get("cell_type") for r in module.rules), [r.id for r in module.rules if not r.when]
    by_source = {}
    for r in module.rules:
        by_source.setdefault(r.source, []).append(r.when["cell_type"])
    assert by_source == {
        "E1": ["K562"],
        "E2": ["CD8_positive__alpha_beta_T_cell"],
        "E3": [measured.CONTEXT_UNKNOWN],
        "E1_measured": ["K562", "WTC11"],
    }
    assert f"# test: rules == {len(module.rules)}" in text
    assert "AlphaGenome deletion, CD8-positive, alpha-beta T cell" in text  # the name itself kept


# --- the split rule of 2026-09-28, applied within each cell -----------------------------------------
def test_a_held_out_cell_is_its_own_marked_rule_and_never_moves_the_training_one():
    r = row_of(pair("AAA", -0.2, "K562"), pair("AAA", -0.8, "WTC11", split=measured.HELDOUT))
    links = sorted(measured.rule_links(r), key=lambda x: x[3])
    assert links == [
        ("AAA", "activates", 0.2, "K562", measured.TRAINING),
        ("AAA", "activates", 0.8, "WTC11", measured.HELDOUT),
    ]
    text, block = program(r)
    rules = {x.when["cell_type"]: x for x in measured_rules(parse(text))}
    assert rules["K562"].strength == 0.2 and not measured.is_heldout(rules["K562"].evidence.source)
    assert measured.is_heldout(rules["WTC11"].evidence.source)
    assert [line for line in block if "targets:" in line] == ["  targets: AAA"]


def test_within_one_cell_the_training_pairs_still_set_the_number():
    r = row_of(pair("AAA", -0.2, "K562"), pair("AAA", -0.9, "K562", split=measured.HELDOUT))
    assert measured.rule_links(r) == [("AAA", "activates", 0.2, "K562", measured.TRAINING)]


def test_the_grammar_says_what_unknown_means_in_a_when_clause():
    text = Path("docs/BIOLANG-GRAMMAR.md").read_text()
    assert "`unknown`" in text and "matches no" in text
