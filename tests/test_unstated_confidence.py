"""A confidence the program did not state is kept apart from a stated 0.0 wherever it travels.

Review item 12 S1 (2026-09-28): the parser gives a missing `confidence:` the value UNSTATED, equal to
0.0 in every calculation, and a stated `confidence: 0.0` exists on purpose (bio.std methylation), so
the value alone cannot carry the difference. These tests follow one small program through every
engine path the census in docs/BIOIR-v0.1.md names: the parser, BioIR JSON (`to_dict` / `from_dict`
through `json`), the BioLang text the `bio` toolchain writes back, and the summaries.

Engine only: this file imports the language and the IR and nothing of the application, so it travels
with the engine package. The application's writers (the `genomeos` CLI, the web endpoints, the forge)
are tested in tests/test_evidence.py and tests/test_unstated_confidence_writers.py.
"""

from __future__ import annotations

import copy
import json
import pickle

from genomeos.ir import (
    RECORDS_UNSTATED,
    UNSTATED,
    UNSTATED_JSON,
    Module,
    confidence_from_json,
    confidence_stated,
    confidence_to_json,
)
from genomeos.lang import parse

PROGRAM = """module test.unstated_paths

param judged_zero = 0.1 {
  evidence: inferred "a guess somebody judged worthless"
  confidence: 0.0
}

param judged_strong = 0.3 { evidence: curated "BioModels"; confidence: 0.9 }

param left_out = 0.2 { evidence: inferred "order of magnitude" }

gene NANOG {
  symbol: NANOG
  evidence: curated "GENCODE v50"
  confidence: 0.8
  transcript NANOG-201 { evidence: curated "GENCODE v50"; confidence: 0.0 }
}

gene SOX2 { symbol: SOX2; evidence: curated "GENCODE v50" }

element E1 {
  class: enhancer
  locus: chr12:7780000-7781000
  targets: NANOG
  evidence: predicted "AlphaGenome deletion" effect -0.2 log2 fold change, probability unavailable
}

rule SOX2 activates NANOG {
  strength: 0.2
  evidence: predicted "AlphaGenome deletion" effect -0.2 log2 fold change, probability unavailable
}

rule NANOG activates SOX2 { strength: 0.3; evidence: curated "a paper"; confidence: 0.0 }
"""

# what the program states, block by block: True = a stated confidence, False = left out
STATED = {
    "judged_zero": True,
    "judged_strong": True,
    "left_out": False,
    "NANOG": True,
    "NANOG-201": True,
    "SOX2": False,
    "E1": False,
    "SOX2 activates NANOG": False,
    "NANOG activates SOX2": True,
}


def _stated(m: Module) -> dict[str, bool]:
    out = {name: confidence_stated(p.confidence) for name, p in m.parameters.items()}
    out |= {e.id: confidence_stated(e.confidence) for e in m.entities.values()}
    out |= {r.id: confidence_stated(r.confidence) for r in m.rules}
    return out


def _values(m: Module) -> dict[str, float]:
    out = {name: float(p.confidence) for name, p in m.parameters.items()}
    out |= {e.id: float(e.confidence) for e in m.entities.values()}
    out |= {r.id: float(r.confidence) for r in m.rules}
    return out


def _json(m: Module) -> Module:
    return Module.from_dict(json.loads(json.dumps(m.to_dict())))


def test_the_parser_keeps_a_stated_zero_and_marks_a_missing_confidence():
    m = parse(PROGRAM)
    got = _stated(m)
    assert {k: got[k] for k in STATED if k != "NANOG-201"} == {
        k: v for k, v in STATED.items() if k != "NANOG-201"
    }
    assert m.parameters["judged_zero"].confidence == 0.0 and m.parameters["left_out"].confidence == 0.0


def test_a_transcript_that_states_zero_keeps_it_rather_than_its_genes():
    """`cconf or conf` treated a stated 0.0 as absent and gave the transcript its gene's 0.8."""
    tx = parse(PROGRAM).entities["NANOG-201"]
    assert confidence_stated(tx.confidence) and tx.confidence == 0.0


def test_bioir_json_keeps_stated_and_unstated_apart():
    m = parse(PROGRAM)
    back = _json(m)
    assert _stated(back) == _stated(m)
    assert _values(back) == _values(m)


def test_bioir_json_writes_an_unstated_confidence_as_null():
    d = json.loads(json.dumps(parse(PROGRAM).to_dict()))
    params = {p["name"]: p for p in d["parameters"]}
    assert params["left_out"]["confidence"] is None
    assert params["judged_zero"]["confidence"] == 0.0


def test_an_enhancer_targets_confidence_survives_json():
    back = _json(parse(PROGRAM))
    assert not confidence_stated(back.entities["E1"].targets[0]["confidence"])


def test_an_ir_object_built_without_a_confidence_is_unstated():
    from genomeos.ir import Gene, Parameter, Rule

    assert not confidence_stated(Gene(id="X", kind="gene").confidence)
    assert not confidence_stated(Rule(id="r", source="X", action="activates", target="X").confidence)
    assert not confidence_stated(Parameter(name="p", value=1.0).confidence)


def test_copy_and_pickle_keep_the_mark():
    for f in (copy.copy, copy.deepcopy, lambda x: pickle.loads(pickle.dumps(x))):
        assert not confidence_stated(f(UNSTATED)) and f(UNSTATED) == 0.0
    m = parse(PROGRAM)
    assert _stated(copy.deepcopy(m)) == _stated(m)


def test_the_confidence_summary_averages_stated_values_only():
    """`confidence_report` pooled an unstated confidence into the mean as a 0.0."""
    m = parse(PROGRAM)
    rep = m.confidence_report()
    assert rep["rule"] == 0.0  # the one stated rule states 0.0; the unstated one is not in the mean
    assert rep["gene"] == 0.8  # SOX2 states none
    assert rep["regulatory_element"] is None  # nothing states one


def test_the_check_report_does_not_list_an_unstated_rule_as_weak():
    from genomeos.lang.tools import check_module

    text = check_module(parse(PROGRAM))
    assert "rules with confidence < 0.5: 1" in text
    assert "SOX2 activates NANOG" not in text.split("rules with confidence < 0.5")[1].split("\n  ")[0]


def test_the_json_form_is_null_and_its_inverse_is_unstated():
    assert UNSTATED_JSON is None and confidence_to_json(UNSTATED) is None
    assert confidence_to_json(0.0) == 0.0 and confidence_stated(confidence_to_json(0.0))
    assert not confidence_stated(confidence_from_json(None)) and confidence_from_json(0.0) == 0.0
    assert confidence_stated(confidence_from_json(0.0))


def test_a_file_written_before_the_mark_reads_its_zero_as_the_stated_zero_it_carries():
    """A BioIR file without `records_unstated_confidence` predates the difference: its 0.0 was written
    the same way whether stated or not. from_dict reads what the file says and does not guess."""
    d = json.loads(json.dumps(parse(PROGRAM).to_dict()))
    assert d.pop(RECORDS_UNSTATED) is True
    for p in d["parameters"]:
        p["confidence"] = 0.0 if p["confidence"] is None else p["confidence"]  # the old writer
    old = Module.from_dict(d)
    assert confidence_stated(old.parameters["left_out"].confidence)  # cannot be told apart: read as written


def test_bio_compile_writes_null_and_the_mark():
    from genomeos.lang.tools import compile_module

    d = json.loads(compile_module(parse(PROGRAM)))
    assert d[RECORDS_UNSTATED] is True
    rules = {r["id"]: r["confidence"] for r in d["rules"]}
    assert rules == {"SOX2 activates NANOG": None, "NANOG activates SOX2": 0.0}
    targets = next(e for e in d["entities"] if e["id"] == "E1")["targets"]
    assert targets[0]["confidence"] is None  # the element's copy, written the same way


def test_the_uncertainty_report_counts_unstated_items_apart():
    from genomeos.runtime.uncertainty import UncertaintyReport, report_for_network

    m = parse(PROGRAM)
    rep = report_for_network(m, m.rules).to_dict()["molecular"]
    # stated: rule NANOG activates SOX2 (0.0 x curated 0.9) and two parameters, inferred 0.0 and curated 0.9
    assert rep["items"] == 5 and rep["unstated"] == 2
    assert rep["confidence"] == round((0.0 + 0.0 * 0.4 + 0.9 * 0.9) / 3, 3)
    only = UncertaintyReport()
    only.add("tissue", UNSTATED, m.rules[0].evidence.kind, "nobody judged this")
    assert only.levels["tissue"].confidence is None and only.levels["tissue"].label == "UNKNOWN"
    assert "none states a confidence" in only.format()


def test_a_trace_line_prints_unstated_rather_than_zero():
    from genomeos.runtime.debugger import TraceLine

    assert TraceLine(1.0, "X", "m", "predicted: m", UNSTATED).format().endswith("conf=unstated")
    assert TraceLine(1.0, "X", "m", "curated: p", 0.0).format().endswith("conf=0.00")
    assert "conf" not in TraceLine(1.0, "X", "m", "runtime default", None).format()


def test_the_bio_confidence_test_is_a_mean_over_stated_rules():
    from genomeos.bio import evaluate

    got = evaluate(parse(PROGRAM), [("confidence", "", ">=", 0.0)])[0]["got"]
    stated = [r.confidence for r in parse(PROGRAM).rules if confidence_stated(r.confidence)]
    assert got == sum(stated) / len(stated)
    none = 'module t\ngene G { evidence: curated "x" }\nrule G activates G { evidence: predicted "m" }\n'
    res = evaluate(parse(none), [("confidence", "", ">=", 0.0)])[0]
    assert res["got"] is None and res["ok"] is False  # nothing stated: the claim cannot be checked


def test_parse_json_and_back_to_biolang_keep_stated_and_unstated():
    """parse -> to_dict -> json -> from_dict -> BioLang text -> parse. No Module-to-BioLang writer
    exists (census, docs/BIOIR-v0.1.md), so the text is written here the way every BioLang writer
    must: the key only when stated, decided by `confidence_to_json`."""
    first = parse(PROGRAM)
    back = _json(first)
    text = ["module test.again"]
    for p in back.parameters.values():
        c = confidence_to_json(p.confidence)
        key = "" if c is None else f"; confidence: {c}"
        text.append(
            f'param {p.name} = {p.value} {{ evidence: {p.evidence.kind.value} "{p.evidence.source}"{key} }}'
        )
    again = parse("\n".join(text))
    assert {k: confidence_stated(v.confidence) for k, v in again.parameters.items()} == {
        k: confidence_stated(v.confidence) for k, v in first.parameters.items()
    }
    assert {k: float(v.confidence) for k, v in again.parameters.items()} == {
        k: float(v.confidence) for k, v in first.parameters.items()
    }
