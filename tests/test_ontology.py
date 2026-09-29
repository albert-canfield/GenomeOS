# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review R7 (2026-09-28): origin, molecular role, activity, target relation and evidence status kept apart.

Before this, `attribution/compile.py` wrote `class: enhancer` on every element, including the 156,925
whose rule inhibits its gene, and each UNKNOWN region carried one role string joining its tier, a
function guess and a reading of constraint. The review's three acceptance tests are here as written,
plus the two constraint rules the pre-registration (docs/ATTRIBUTION.md, R7) adds.
"""

from __future__ import annotations

import gzip
import json

import pytest

from genomeos.attribution.compile import REPRESSION_ALTERNATIVES, compile_chromosome
from genomeos.attribution.measured import Layer
from genomeos.ir.model import Module
from genomeos.lang import grammar
from genomeos.lang.parser import BioLangError, parse

GROUP = list(REPRESSION_ALTERNATIVES)


def _fixture(tmp_path):
    blocks = [
        {  # a regulatory block over a LINE: repeat-derived and regulatory at once
            "start": 0,
            "end": 1000,
            "length": 1000,
            "class": "regulatory",
            "class_evidence": "curated",
            "phylop": {"bases": 1000, "fraction_above": 0.2},
            "guess": {"tier": "regulatory", "label": "regulatory elements", "confidence": 0.7},
        },
        {  # an unconstrained fossil: no sign of selection, which is not a sign of no function
            "start": 5000,
            "end": 6000,
            "length": 1000,
            "class": "interspersed_repeat_SINE",
            "class_evidence": "curated",
            "phylop": {"bases": 1000, "fraction_above": 0.0},
            "guess": {
                "tier": "fossil",
                "label": "transposable-element fossil, unconstrained",
                "confidence": 0.8,
            },
        },
        {  # unconstrained unique sequence, the budget's "best guess neutral"
            "start": 7000,
            "end": 8000,
            "length": 1000,
            "class": "unique_intergenic",
            "class_evidence": "inferred",
            "phylop": {"bases": 1000, "fraction_above": 0.01},
            "guess": {"tier": "neutral", "label": "unconstrained unique sequence", "confidence": 0.6},
        },
    ]
    (tmp_path / "budget_chrT.json").write_text(
        json.dumps({"chrom": "chrT", "unknown_bp": 3000, "constrained_fraction": 0.1, "blocks": blocks})
    )
    (tmp_path / "unknown_chrT.json").write_text(
        json.dumps(
            {
                "blocks": [
                    {"start": 0, "features": {"ccre": {"dELS": 3, "CTCF-only": 1}}},
                    {"start": 7000, "features": {"interspersed_coverage": 0.0, "repeat_coverage": {}}},
                ]
            }
        )
    )
    (tmp_path / "domains_chrT.json").write_text(
        json.dumps({"domains": [{"id": "chrT:D1", "start": 0, "end": 100000}]})
    )
    with gzip.open(tmp_path / "ccres_chrT.bed.gz", "wt") as fh:
        fh.write("# fixture\n")
        fh.write("chrT\t100\t300\tE_line\tdELS\t0\n")
        fh.write("chrT\t2000\t2200\tE_up\tdELS\t1\n")
    with gzip.open(tmp_path / "rmsk_chrT.bed.gz", "wt") as fh:
        fh.write("50\t900\tLINE\tL2\tL2a\t0.200\n")
        fh.write("5000\t6000\tSINE\tAlu\tAluY\t0.050\n")

    def element(eid, start, end, lfc, action):
        return {
            "id": eid,
            "start": start,
            "end": end,
            "domain": "chrT:D1",
            "verdict_coding": "agrees with nearest TSS in domain",
            "predicted_coding": {
                "gene": "G1",
                "action": action,
                "log2_fold_change": lfc,
                "tissue": "K562",
                "strength": "weak",
                "confidence": abs(lfc),
            },
        }

    els = [element("E_line", 100, 300, -0.6, "activates"), element("E_up", 2000, 2200, 0.5, "represses")]
    (tmp_path / "enhancer_targets_chrT.json").write_text(json.dumps({"elements": els}))


@pytest.fixture
def module(tmp_path):
    _fixture(tmp_path)
    text = compile_chromosome("chrT", tmp_path, layer=Layer(chrom="chrT"))
    return text, parse(text, "chrT")


def axes(m, eid):
    return m.entities[eid].attrs["ontology"]


def test_a_repeat_derived_regulatory_element_keeps_both(module):
    _, m = module
    a = axes(m, "E_line")
    assert a["origin"] == [["repeat_derived/LINE"]]
    assert ["enhancer_like"] in a["molecular_role"]
    assert m.entities["E_line"].cls == "enhancer"
    region = axes(m, "U_chrT_0")
    assert region["origin"] == [["repeat_derived/LINE"]]
    assert region["molecular_role"] == [["enhancer_like"], ["insulator_like"]]
    assert "registry_biochemical" in [g[0] for g in region["evidence_status"]]


def test_a_predicted_increase_does_not_force_a_silencer_label(module):
    _, m = module
    el = m.entities["E_up"]
    a = el.attrs["ontology"]
    assert el.cls == "enhancer"  # the registry summary, not the direction
    assert a["activity"] == [["represses_target"]]
    assert ["silencer"] not in a["molecular_role"]
    assert GROUP in a["molecular_role"]  # left open, with unknown among the alternatives
    assert ["insulator_like"] in a["molecular_role"]  # CTCF-bound: overlapping roles both hold
    rule = next(r for r in m.rules if r.source == "E_up")
    assert rule.action.value == "inhibits"


def test_unresolved_alternatives_survive_serialisation(module):
    _, m = module
    back = Module.from_dict(json.loads(json.dumps(m.to_dict())))
    for eid in ("E_up", "E_line", "U_chrT_0", "U_chrT_5000", "U_chrT_7000"):
        assert back.entities[eid].attrs["ontology"] == m.entities[eid].attrs["ontology"]
    assert GROUP in back.entities["E_up"].attrs["ontology"]["molecular_role"]


def test_constraint_is_selection_and_its_absence_is_not_no_function(module):
    _, m = module
    for eid in ("U_chrT_5000", "U_chrT_7000"):
        a = axes(m, eid)
        assert a["molecular_role"] == [["unknown"]]
        assert ["selection_not_detected"] in a["evidence_status"]
    assert ["under_selection"] in axes(m, "U_chrT_0")["evidence_status"]
    assert axes(m, "U_chrT_5000")["origin"] == [["repeat_derived/SINE"]]
    assert axes(m, "U_chrT_7000")["origin"] == [["unique"]]
    # no value in the vocabulary means "no function"
    assert not {v for v in grammar.AXES["molecular_role"] if "neutral" in v or "nonfunction" in v}
    # the predicted elements carry no constraint, so they say so rather than guess
    assert ["selection_not_measured"] in axes(m, "E_line")["evidence_status"]


def test_the_counts_the_program_tests_do_not_move(module):
    text, m = module
    assert "# test: unknowns == 0" in text and len(m.rules) == 2


def test_a_value_outside_the_vocabulary_is_refused():
    with pytest.raises(BioLangError, match="molecular_role has no value 'silencr'"):
        parse("element X { class: enhancer; molecular_role: enhancer_like, silencr|unknown }", "t")
    m = parse(
        "element X { origin: repeat_derived/LTR; molecular_role: enhancer_like, silencer|unknown }", "t"
    )
    assert m.entities["X"].attrs["ontology"] == {
        "origin": [["repeat_derived/LTR"]],
        "molecular_role": [["enhancer_like"], ["silencer", "unknown"]],
    }


def test_a_block_without_axes_parses_as_before():
    m = parse("element X { class: enhancer }\nregion R { role: unknown }", "t")
    assert "ontology" not in m.entities["X"].attrs and "ontology" not in m.entities["R"].attrs
