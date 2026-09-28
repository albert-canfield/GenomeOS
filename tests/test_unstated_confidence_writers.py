# SPDX-License-Identifier: AGPL-3.0-or-later
"""The application's writers keep an unstated confidence apart from a stated 0.0 (item 12 S1).

The census in docs/BIOIR-v0.1.md (2026-09-28) names every path an IR confidence travels. The engine's
paths are tested in tests/test_unstated_confidence.py; these are the application's: the `genomeos
compile` and `genomeos check` commands, the web endpoints that send an IR confidence as JSON (compile,
cell, debug), an organism experiment's result, and the forge's experiment.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.ir import Module, confidence_stated
from genomeos.lang import parse, parse_file

CENSUS = pytest.mark.xfail(
    strict=True,
    reason="census 2026-09-28 (item 12 S1): this writer loses an unstated confidence; the build fixes it",
)

ROOT = Path(__file__).resolve().parent.parent

PROGRAM = """module test.writers

param mrna_half_life = 2 { evidence: inferred "order of magnitude" }

param judged_zero = 0.1 { evidence: inferred "judged worthless"; confidence: 0.0 }

gene A { basal: 1; max: 5; produces: Ap; evidence: curated "GENCODE v50"; confidence: 0.9 }

gene B { basal: 0; max: 5; produces: Bp; evidence: curated "GENCODE v50" }

protein Ap { evidence: curated "UniProt"; confidence: 0.9 }

protein Bp { evidence: curated "UniProt" }

rule Ap activates B { strength: 0.5; threshold: 1; hill: 2; evidence: predicted "a model" }

rule Bp inhibits A { strength: 0.2; threshold: 1; hill: 2; evidence: curated "a paper"; confidence: 0.0 }

cell_type T { expresses: A, B; evidence: curated "an atlas" }
"""


@CENSUS
def test_genomeos_compile_writes_null_for_unstated_and_check_reads_it_back(tmp_path, capsys):
    from genomeos.cli import main

    src = tmp_path / "w.bio"
    src.write_text(PROGRAM)
    out = tmp_path / "w.json"
    assert main(["compile", str(src), "-o", str(out)]) == 0
    d = json.loads(out.read_text())
    assert d["records_unstated_confidence"] is True
    params = {p["name"]: p for p in d["parameters"]}
    assert params["mrna_half_life"]["confidence"] is None and params["judged_zero"]["confidence"] == 0.0
    rules = {r["id"]: r["confidence"] for r in d["rules"]}
    assert rules == {"A->Ap": 0.9, "B->Bp": None, "Ap activates B": None, "Bp inhibits A": 0.0}
    back = {r.id: confidence_stated(r.confidence) for r in Module.from_dict(d).rules}
    assert back == {"A->Ap": True, "B->Bp": False, "Ap activates B": False, "Bp inhibits A": True}
    capsys.readouterr()
    assert main(["check", str(out)]) == 0
    from_json = capsys.readouterr().out
    assert main(["check", str(src)]) == 0
    assert capsys.readouterr().out == from_json  # the JSON says what the program says


@CENSUS
def test_genomeos_check_lists_only_stated_weak_rules(tmp_path, capsys):
    from genomeos.cli import main

    src = tmp_path / "w.bio"
    src.write_text(PROGRAM)
    assert main(["check", str(src)]) == 0
    out = capsys.readouterr().out
    assert "rules with confidence < 0.5: 1\n    Bp inhibits A" in out
    assert "rules that state no confidence: 2 (not counted as weak)" in out
    assert "gene       ██████████████████   0.90  (1 of 2 unstated)" in out
    assert "cell_type" in out and "none stated (1 unstated)" in out


@CENSUS
def test_the_web_compile_report_keeps_unstated_apart():
    from genomeos.web.server import Api

    r = Api(ROOT).compile(PROGRAM)
    assert r["ok"]
    rep = r["report"]
    assert rep["confidence"]["gene"] == 0.9 and rep["confidence"]["cell_type"] is None
    assert rep["confidence"]["rule"] == 0.45  # (0.9 + 0.0) / 2: the two unstated rules are not zeros
    assert rep["confidence_counts"]["rule"] == {"stated": 2, "unstated": 2}
    assert [w["id"] for w in rep["weak_rules"]] == ["Bp inhibits A"] and rep["unstated_rules"] == 2
    module = json.loads(json.dumps(r["module"]))
    got = {x["id"]: x["confidence"] for x in module["rules"]}
    assert got["Ap activates B"] is None and got["Bp inhibits A"] == 0.0


@CENSUS
def test_the_web_cell_view_sends_null_for_an_unstated_confidence(tmp_path):
    from genomeos.web.server import Api

    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "w.bio").write_text(PROGRAM)
    out = json.loads(json.dumps(Api(tmp_path).cell("data/w.bio", "T")))
    assert out["cell_type"]["confidence"] is None
    rules = {r["id"]: r["confidence"] for r in out["rules"]}
    assert rules["Ap activates B"] is None and rules["Bp inhibits A"] == 0.0 and rules["A->Ap"] == 0.9
    edges = {(e["source"], e["target"]): e["confidence"] for e in out["graph"]["edges"]}
    assert edges[("Ap", "B")] is None and edges[("Bp", "A")] == 0.0


@CENSUS
def test_the_web_debugger_sends_null_for_an_unstated_term():
    from genomeos.web.server import Api

    r = Api(ROOT).debug(PROGRAM, [], {"A": 1.0}, until=1.0, explain=["B.mRNA"])
    lines = json.loads(json.dumps(r))["explain"]["B.mRNA"]
    by = {x["message"].split("  ", 1)[1]: x["confidence"] for x in lines}
    assert by["basal transcription of B"] is None  # a runtime default: no evidence, no confidence
    assert next(v for k, v in by.items() if k.startswith("Ap activates B")) is None
    assert by["mRNA decay"] is None  # the program's mrna_half_life, which states none


@CENSUS
def test_an_experiment_result_writes_null_for_an_experiment_that_states_none():
    from genomeos.organism.experiment import run_experiments

    src = """
module toy.ko
import bio.std.development
organism T { root: Z; cell_type: Zygote; factors: F }
cell_type WithF { parent: Blastomere }
cell_type WithoutF { parent: Blastomere }
decision z { action: divide; when: cell = Z; daughters: A, B; after: 1 min; asymmetric: F -> A }
decision with { action: differentiate; when: cell = A, F = present; to: WithF }
decision without { action: differentiate; when: cell = A, F = absent; to: WithoutF }
experiment unjudged { knockout: F; until: 10 min }
experiment judged { knockout: F; until: 10 min; evidence: inferred "a guess"; confidence: 0.0 }
"""
    unjudged, judged = run_experiments(parse(src))
    assert json.loads(json.dumps(unjudged.to_dict()))["evidence"]["confidence"] is None
    assert json.loads(json.dumps(judged.to_dict()))["evidence"]["confidence"] == 0.0


@CENSUS
def test_the_forge_experiment_agrees_with_the_text_it_emits():
    """BioForge states no confidence on its experiment; the parsed text and the returned object used
    to disagree, because the object carried a plain 0.0 that `confidence_stated` read as judged."""
    from genomeos.organism.forge import UNSTATED_CONFIDENCE, run_designs

    assert not confidence_stated(UNSTATED_CONFIDENCE) and UNSTATED_CONFIDENCE == 0.0
    results = run_designs(parse_file(str(ROOT / "data/organisms/celegans/designs.bio")), None, seed=0)
    r = next(x for x in results if x.to_experiment() is not None)
    ex = r.to_experiment()
    reparsed = parse(f"module t\n{r.to_bio()}").experiments[0]
    assert confidence_stated(ex.confidence) is confidence_stated(reparsed.confidence) is False
