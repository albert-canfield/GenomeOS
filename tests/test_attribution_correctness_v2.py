# SPDX-License-Identifier: AGPL-3.0-or-later
"""The direction rule, versioned (registered 2026-09-29, lane-judge2): v1 is S4's rule and reproduces its
committed result; v2 judges a direction only in the stated cell and keeps every other cell as a finding."""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import itertools
import json
import re
import sys
from argparse import Namespace
from pathlib import Path

import pytest

from genomeos.attribution import correctness as co
from genomeos.attribution import holdout as ho
from genomeos.attribution import measured as ms

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location("s4_correctness_run", ROOT / "scripts/s4_correctness_run.py")
s4 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s4)

V1_RESULT = ROOT / "data/results/attribution_correctness.json"
V2_RESULT = ROOT / "data/results/attribution_correctness_v2.json"
E = ("chr1", 10_000, 10_400)
CELLS = ("K562", "HepG2", "placenta", "", "unknown")
KINDS = ("crispri_decrease", "crispri_increase", "crispri_null_well_powered", "crispri_null_underpowered")


def claim(axis: str, value: str, gene: str = "G1", cell: str = "K562") -> co.Claim:
    return co.Claim("EH1", *E, axis, value, gene, cell)


def obs(kind: str, cell: str = "K562", gene: str = "G1") -> co.Observation:
    return co.Observation(kind, f"crispri:{cell or 'x'}", gene, cell, "training")


def both(c: co.Claim, seen: list[co.Observation]) -> tuple[co.Verdict, co.Verdict]:
    return co.verdict_of(c, seen, co.RULE_V1), co.verdict_of(c, seen, co.RULE_V2)


# --- the rule ----------------------------------------------------------------------------------------
def test_v2_stays_selectable_and_an_unknown_rule_is_refused():
    # v2 was the default from its registration until v3's (tests/test_attribution_correctness_v3.py)
    assert co.RULES[:2] == (co.RULE_V1, co.RULE_V2) and co.DEFAULT_RULE != co.RULE_V2
    c = claim(co.ACTIVITY, "activates_target", cell="placenta")
    assert co.verdict_of(c, [obs("crispri_decrease")], co.RULE_V2).rule == co.RULE_V2
    with pytest.raises(ValueError):
        co.verdict_of(c, [], rule="v4")
    with pytest.raises(ValueError):
        co.judge([], ho.Labels("x", lambda u, e: None, frozenset()), units={}, sources=[], rule="v4")


def test_agreement_in_another_cell_no_longer_establishes_and_disagreement_no_longer_refutes():
    fell_in_k562 = [obs("crispri_decrease", "K562")]
    agree = claim(co.ACTIVITY, "activates_target", cell="placenta")
    disagree = claim(co.ACTIVITY, "represses_target", cell="placenta")
    a1, a2 = both(agree, fell_in_k562)
    d1, d2 = both(disagree, fell_in_k562)
    assert (a1.verdict, d1.verdict) == (co.CORRECT, co.INCORRECT)  # v1, S4 as registered
    for v, finding in ((a2, co.AGREES), (d2, co.DISAGREES)):
        assert (v.verdict, v.reason, v.detail) == (
            co.NOT_JUDGED,
            co.NOT_ASSESSED,
            "the gene responded only in other cells",
        )
        assert v.deciding == () and v.refutable is False
        assert v.cross_cell == ((fell_in_k562[0], finding),)


def test_a_claim_with_no_stated_cell_is_not_assessed_in_this_context():
    for cell in ("", "unknown"):
        v = co.verdict_of(claim(co.ACTIVITY, "activates_target", cell=cell), [obs("crispri_decrease")])
        assert (v.verdict, v.reason, v.detail) == (co.NOT_JUDGED, co.NOT_ASSESSED, "the claim states no cell")
        assert [f for _, f in v.cross_cell] == [co.AGREES]


def test_the_stated_cell_decides_and_other_cells_stay_beside():
    seen = [obs("crispri_decrease", "K562"), obs("crispri_increase", "HepG2")]
    v = co.verdict_of(claim(co.ACTIVITY, "activates_target", cell="K562"), seen)
    assert (v.verdict, v.deciding, v.refutable) == (co.CORRECT, (seen[0],), True)
    assert v.cross_cell == ((seen[1], co.DISAGREES),)
    both_here = [obs("crispri_decrease", "K562"), obs("crispri_increase", "K562")]
    v = co.verdict_of(claim(co.ACTIVITY, "activates_target", cell="K562"), both_here)
    assert v.verdict == co.MODEL_INADEQUATE and v.refutable is True and v.cross_cell == ()
    # v1 called a split across two other cells unresolved; v2 does not assess it and keeps both
    v = co.verdict_of(claim(co.ACTIVITY, "activates_target", cell="WTC11"), seen)
    assert v.reason == co.NOT_ASSESSED
    assert sorted(f for _, f in v.cross_cell) == [co.AGREES, co.DISAGREES]


@pytest.mark.parametrize("value", co.DIRECTION)
@pytest.mark.parametrize("cell", CELLS)
def test_v2_never_decides_a_direction_from_another_cell(value, cell):
    """Every combination of CRISPRi outcomes in two cells: v2 decides only with stated-cell observations;
    a v1 verdict decided in another cell only is not assessed under v2; everything else is v1's."""
    outcomes = [(k, c) for k in KINDS for c in ("K562", "HepG2")]
    c = claim(co.ACTIVITY, value, cell=cell)
    for n in range(3):
        for combo in itertools.combinations(outcomes, n):
            seen = [obs(k, oc) for k, oc in combo]
            v1, v2 = both(c, seen)
            assert v2.rule == co.RULE_V2 and v2.cross_cell is not None
            assert all(co.stated(cell) and co.norm_cell(o.cell) == co.norm_cell(cell) for o in v2.deciding)
            assert v2.refutable == (v2.verdict in co.JUDGED)
            here = [o for o in v1.deciding if co.stated(cell) and co.norm_cell(o.cell) == co.norm_cell(cell)]
            if v1.verdict in co.JUDGED and not here:  # decided only in another cell
                assert (v2.verdict, v2.reason) == (co.NOT_JUDGED, co.NOT_ASSESSED), combo
                assert {o for o, _ in v2.cross_cell} >= set(v1.deciding)
            elif v1.verdict in co.JUDGED and len(here) == len(v1.deciding):  # decided in the stated cell
                assert (v2.verdict, v2.deciding) == (v1.verdict, v1.deciding), combo
            elif v1.verdict in co.JUDGED:  # v1 mixed the stated cell with others (unresolved)
                assert set(v2.deciding) <= set(here), combo
            else:
                assert (v2.verdict, v2.reason, v2.detail) == (v1.verdict, v1.reason, v1.detail), combo


def test_nothing_but_the_direction_moves_between_the_rules():
    seen = [
        obs("crispri_decrease", "K562"),
        obs("crispri_null_well_powered", "HepG2"),
        co.Observation("lentimpra_active", "lentimpra:K562", cell="K562"),
        co.Observation("vista_positive", "vista", cell=co.VISTA_CONTEXT),
        co.Observation("registry_biochemical", "encode", value="enhancer_like"),
    ]
    others = [
        claim(co.TARGET, "G1", cell="HepG2"),
        claim(co.TARGET, "G1", cell="placenta"),
        claim(co.CONTEXT, "HepG2", cell="HepG2"),
        claim(co.CONTEXT, "placenta", cell="placenta"),
        claim(co.ACTIVITY, "active_in_reporter", gene="", cell="HepG2"),
        claim(co.ACTIVITY, "inactive_in_reporter", gene="", cell=""),
        claim(co.ACTIVITY, "no_effect_measured", cell="HepG2"),
        claim(co.ROLE, "enhancer_like", gene="", cell=""),
        claim(co.ORIGIN, "unique", gene="", cell=""),
    ]
    for c in others:
        v1, v2 = both(c, seen)
        assert s4._pair(v1) == s4._pair(v2), c
        assert v2.cross_cell is None
        if c.axis in (co.TARGET, co.CONTEXT):
            assert v2.refutable == v1.refutable
        else:  # not computed: written as null under v2, never false
            assert v1.refutable is False and v2.refutable is None


def test_a_v1_verdict_and_report_are_written_exactly_as_s4_wrote_them():
    c = claim(co.ACTIVITY, "activates_target", cell="placenta")
    d = co.verdict_of(c, [obs("crispri_decrease")], co.RULE_V1).to_dict()
    assert "rule" not in d and "cross_cell" not in d and d["refutable"] is False
    t = co.AxisTally()
    assert (
        co.NOT_ASSESSED not in t.to_dict()["not_judged_by_reason"]
        and "cross_cell_findings" not in t.to_dict()
    )
    assert "refutable" not in t.accuracy(co.ACTIVITY).beside
    assert "refutable" in co.AxisTally(rule=co.RULE_V2).accuracy(co.ACTIVITY).beside


def test_the_cross_cell_finding_survives_serialisation():
    c = claim(co.ACTIVITY, "represses_target", cell="Whole_Blood")
    v = co.verdict_of(c, [obs("crispri_decrease", "K562")], co.RULE_V2)
    d = json.loads(json.dumps(v.to_dict()))
    assert d["rule"] == co.RULE_V2 and d["reason"] == co.NOT_ASSESSED
    assert d["cross_cell"] == [{**obs("crispri_decrease", "K562").to_dict(), "finding": co.DISAGREES}]
    assert co.Verdict.from_dict(d) == v
    role = co.verdict_of(claim(co.ROLE, "enhancer_like", gene="", cell=""), [])
    assert co.Verdict.from_dict(json.loads(json.dumps(role.to_dict()))) == role


def test_judge_keeps_the_not_assessed_claims_and_their_findings(tmp_path):
    units = {
        "crispri:A": (
            ho.Unit(
                source="crispri:A",
                chrom=E[0],
                start=E[1],
                end=E[2],
                outcome=ms.DECREASE,
                gene="G1",
                cell="K562",
            ),
        )
    }
    claims = [
        claim(co.ACTIVITY, "represses_target", cell="Whole_Blood"),
        claim(co.TARGET, "G1", cell="Whole_Blood"),
    ]
    labels = ho.Labels("synthetic", lambda u, e: None, frozenset({ho.MODEL}))
    r1 = co.judge(claims, labels, units=units, sources=["crispri:A"], rule=co.RULE_V1)
    r2 = co.judge(claims, labels, units=units, sources=["crispri:A"], rule=co.RULE_V2)
    assert [v.verdict for v in r1.judged] == [co.INCORRECT, co.CORRECT] and "rule" not in r1.to_dict()
    assert [v.claim.axis for v in r2.judged] == [co.TARGET]
    assert [v.reason for v in r2.not_assessed] == [co.NOT_ASSESSED]
    act = r2.axes[co.ACTIVITY].to_dict()
    assert act["not_judged_by_reason"][co.NOT_ASSESSED] == 1
    assert act["cross_cell_findings"] == {co.NOT_ASSESSED: {co.DISAGREES: 1}}
    d = json.loads(json.dumps(r2.to_dict()))
    assert d["rule"] == co.RULE_V2 and len(d["not_assessed"]) == 1
    back = co.Report.from_dict(d)
    assert back.to_dict() == d and back.rule == co.RULE_V2
    assert co.Report.from_dict(json.loads(json.dumps(r1.to_dict()))).to_dict() == r1.to_dict()


# --- every historical caller pins v1 ----------------------------------------------------------------
PINNED = ("scripts/prior_only_test.py", "scripts/pilot_biological_gate.py", "scripts/repression_trace.py")


@pytest.mark.parametrize("path", PINNED)
def test_the_scripts_that_reproduce_v1_results_pin_v1(path):
    calls = [
        n
        for n in ast.walk(ast.parse((ROOT / path).read_text()))
        if isinstance(n, ast.Call)
        and isinstance(n.func, ast.Attribute)
        and n.func.attr in ("judge", "verdict_of")
    ]
    assert calls, path
    for n in calls:
        args = [*n.args, *(k.value for k in n.keywords if k.arg == "rule")]
        text = " ".join(ast.unparse(a) for a in args)
        assert "RULE_V1" in text or re.search(r"\bv1\b", text), (path, ast.unparse(n))


def test_the_s4_script_pins_v1_unless_asked_for_v2():
    src = (ROOT / "scripts/s4_correctness_run.py").read_text()
    assert "report = judged_under(co.RULE_V1, a, units, refs)" in src
    assert s4.NAMES == {co.RULE_V1: "attribution_correctness", co.RULE_V2: "attribution_correctness_v2"}


def test_the_document_carries_the_registration_word_for_word():
    doc = (ROOT / "docs/ATTRIBUTION.md").read_text()
    start = doc.index("## The direction rule, versioned")
    section = re.sub(r"\s+", " ", re.sub(r"\n> ?", "\n", doc[start : doc.index("\n## ", start + 5)]))
    for text in (
        co.RULE_DECISION.split(": ", 1)[1],
        co.DIRECTION_RULE[co.RULE_V2],
        co.REFUTABLE_RULE[co.RULE_V2],
        co.RULE_RECORD,
        co.SEEN_BEFORE_V2,
        co.EXPECTED_V2,
    ):
        assert re.sub(r"\s+", " ", text) in section
    assert co.DIRECTION_RULE[co.RULE_V1].endswith(co.AXIS_RULE[co.ACTIVITY])


# --- the committed re-judge --------------------------------------------------------------------------
needs_v2 = pytest.mark.skipif(
    not V2_RESULT.exists(), reason="the v2 re-judge is written once, after the code"
)


@pytest.fixture(scope="module")
def v2() -> dict:
    return json.loads(V2_RESULT.read_text())


@pytest.fixture(scope="module")
def v1() -> dict:
    return json.loads(V1_RESULT.read_text())


@needs_v2
def test_the_v2_result_names_its_rule_and_its_stamp_is_clean(v2):
    m = v2["result_manifest"]
    assert v2["rule"] == co.RULE_V2 and m["parameters"]["judge_rule"] == co.RULE_V2
    assert m["complete"] and m["code"]["dirty"] is False
    assert v2["expected_before_the_run"] == co.EXPECTED_V2 and v2["alphagenome_requests"] == 0


@needs_v2
def test_the_historical_result_is_kept_unchanged_beside_it(v2):
    h = v2["historical_result"]
    assert h["v1_rerun_reproduces_it"] is True and h["rule"] == co.RULE_V1
    read = {i["path"]: i["sha256"] for i in v2["result_manifest"]["inputs"]}
    assert (
        hashlib.sha256(V1_RESULT.read_bytes()).hexdigest()
        == read["data/results/attribution_correctness.json"]
    )


@needs_v2
def test_v1_reproduces_the_committed_s4_counts(v2, v1):
    for axis in co.AXES:
        row, t = v2["side_by_side"][axis][co.RULE_V1], v1["axes"][axis]
        assert {k: row[k] for k in co.VERDICTS} == t["verdicts"]
        assert {r: row["not_judged_by_reason"][r] for r in co.REASONS} == t["not_judged_by_reason"]
        assert row["not_judged_by_reason"][co.NOT_ASSESSED] == 0 and row["claims"] == t["claims"]
    assert v1["axes"][co.ACTIVITY]["verdicts"][co.CORRECT] == 93
    assert v1["axes"][co.ACTIVITY]["verdicts"][co.INCORRECT] == 5


@needs_v2
def test_the_59_move_and_nothing_else(v2):
    sbs = v2["side_by_side"]
    for axis in (co.ORIGIN, co.ROLE, co.TARGET, co.CONTEXT):
        assert sbs[axis][co.RULE_V1] == sbs[axis][co.RULE_V2], axis
    a1, a2 = sbs[co.ACTIVITY][co.RULE_V1], sbs[co.ACTIVITY][co.RULE_V2]
    assert (a2["judged"], a2[co.CORRECT], a2[co.INCORRECT]) == (39, 39, 0)
    assert a2["not_judged_by_reason"][co.NOT_ASSESSED] == 59
    assert a2[co.NOT_JUDGED] - a1[co.NOT_JUDGED] == 59
    for r in co.REASONS:
        assert a2["not_judged_by_reason"][r] == a1["not_judged_by_reason"][r]
    assert (a1["refutable"], a2["refutable"]) == (None, 39)
    m = v2["moved"]
    assert m["claims"] == 59 and m["by_axis"] == {co.ACTIVITY: 59}
    assert m["v1_verdicts"] == {co.CORRECT: 54, co.INCORRECT: 5}
    assert m["v1_decided_in"] == {"another cell only": 59}
    assert m["v2_verdicts"] == {f"{co.NOT_JUDGED}:{co.NOT_ASSESSED}": 59}
    assert m["unexplained_count_changes"] == []


@needs_v2
def test_the_cross_cell_finding_is_kept_on_every_moved_claim(v2):
    said = {}
    for e in v2["moved"]["each"]:
        cross = e[co.RULE_V2]["cross_cell"]
        assert cross and {f["finding"] for f in cross} in ({co.AGREES}, {co.DISAGREES})
        assert {(o["kind"], o["cell"]) for o in e[co.RULE_V1]["deciding"]} <= {
            (o["kind"], o["cell"]) for o in cross
        }
        said[e[co.RULE_V1]["verdict"]] = said.get(e[co.RULE_V1]["verdict"], set()) | {cross[0]["finding"]}
    assert said == {co.CORRECT: {co.AGREES}, co.INCORRECT: {co.DISAGREES}}
    assert v2["moved"]["v2_cross_cell"] == {co.AGREES: 54, co.DISAGREES: 5}
    assert len(v2["not_assessed"]) == 59


@needs_v2
def test_id1_is_not_assessed_and_its_k562_disagreement_is_kept(v2):
    f = v2["separate_findings"]["ID1"]
    (v,) = f["v2_verdicts"]
    assert (v["element"], v["value"], v["cell"]) == ("EH38E3426791", "represses_target", "Whole_Blood")
    assert (v["verdict"], v["reason"]) == (co.NOT_JUDGED, co.NOT_ASSESSED)
    assert [(o["kind"], o["cell"], o["source"], o["finding"]) for o in v["cross_cell"]] == [
        ("crispri_decrease", "K562", "crispri:Gasperini2019", co.DISAGREES)
    ]
    assert f["repression_trace_class_c"].startswith("model on K562: by_cell[K562] +0.1433")
    assert "does not directly refute the compiled whole-blood claim" in f["reading"]


# --- the live judge on the committed claims and sources (local data only) ----------------------------
live = pytest.mark.skipif(
    not (ROOT / ho.COMPILED_DIR).exists() or not (ROOT / "data/knowledge/crispri").exists(),
    reason="the compiled programs and the CRISPRi benchmark are local and untracked",
)


@pytest.fixture(scope="module")
def judged() -> dict:
    with pytest.MonkeyPatch.context() as mp:
        mp.chdir(ROOT)
        units, refs = ho.all_units(), ho.crispri_references()
        a = Namespace(compiled=ho.COMPILED_DIR, chroms=None)
        out = {"units": units, co.RULE_V1: s4.judged_under(co.RULE_V1, a, units, refs)}
        if V2_RESULT.exists():
            out[co.RULE_V2] = s4.judged_under(co.RULE_V2, a, units, refs)
        return out


@live
def test_v1_reproduces_the_committed_s4_result_exactly(judged, v1):
    rerun = json.loads(json.dumps(s4.v1_payload(judged[co.RULE_V1], judged["units"], 0.0), default=str))
    drop = set(s4.RUN_KEYS)
    assert s4._text(rerun, drop) == s4._text(v1, drop)


@live
@needs_v2
def test_v2_reproduces_the_committed_re_judge(judged, v2):
    body = judged[co.RULE_V2].to_dict()
    for k in ("rule", "quantities", "also_reported", "axes", "not_claims", "judged", "not_assessed"):
        assert body[k] == v2[k], k


# --- the v2 run, end to end on a synthetic program ---------------------------------------------------
SYNTHETIC = """\
element EH1 {
  class: enhancer
  locus: chr21:10000-10400
  origin: unique
  molecular_role: enhancer_like
  activity: represses_target
  target_relation: predicted_deletion_target
  evidence_status: predicted_model
}
rule EH1 inhibits G1 { strength: 0.2; when: cell_type = Whole_Blood; evidence: predicted "AlphaGenome" }
element EH2 {
  class: enhancer
  locus: chr21:20000-20400
  origin: unique
  molecular_role: enhancer_like
  activity: activates_target
  target_relation: predicted_deletion_target
  evidence_status: predicted_model
}
rule EH2 activates G2 { strength: 0.2; when: cell_type = K562; evidence: predicted "AlphaGenome" }
"""


def test_the_v2_run_reports_both_rules_and_every_moved_claim(tmp_path, monkeypatch):
    (tmp_path / "noncoding_chr21.bio").write_text(SYNTHETIC)

    def unit(start: int, gene: str) -> ho.Unit:
        return ho.Unit(
            source="crispri:A",
            chrom="chr21",
            start=start,
            end=start + 400,
            outcome=ms.DECREASE,
            gene=gene,
            cell="K562",
        )

    monkeypatch.setattr(s4.ho, "all_units", lambda: {"crispri:A": (unit(10_000, "G1"), unit(20_000, "G2"))})
    monkeypatch.setattr(s4.ho, "crispri_references", lambda: {})
    saved = {}
    monkeypatch.setattr(
        s4, "save_result", lambda name, payload, where, manifest: saved.update(name=name, p=payload)
    )
    monkeypatch.setattr(s4, "manifest", lambda programs, params, also=(): {"parameters": params})
    s4.main(["--rule", "v2", "--compiled", str(tmp_path), "--results", str(tmp_path)])
    p = saved["p"]
    assert saved["name"] == "attribution_correctness_v2" and p["rule"] == co.RULE_V2
    assert p["historical_result"]["v1_rerun_reproduces_it"] is None  # not the registered run: not compared
    act = p["side_by_side"][co.ACTIVITY]
    assert (act[co.RULE_V1][co.INCORRECT], act[co.RULE_V1][co.CORRECT]) == (1, 1)
    assert (act[co.RULE_V2][co.CORRECT], act[co.RULE_V2]["not_judged_by_reason"][co.NOT_ASSESSED]) == (1, 1)
    assert p["side_by_side"][co.TARGET][co.RULE_V1] == p["side_by_side"][co.TARGET][co.RULE_V2]
    m = p["moved"]
    assert (m["claims"], m["v1_decided_in"], m["unexplained_count_changes"]) == (
        1,
        {"another cell only": 1},
        [],
    )
    (e,) = m["each"]
    assert e["claim"]["element"] == "EH1" and e[co.RULE_V2]["cross_cell"][0]["finding"] == co.DISAGREES
    assert [v["element"] for v in p["not_assessed"]] == ["EH1"]


def test_the_run_without_a_rule_writes_s4s_v1_result_as_before(tmp_path, monkeypatch):
    (tmp_path / "noncoding_chr21.bio").write_text(SYNTHETIC)
    unit = ho.Unit(
        source="crispri:A",
        chrom="chr21",
        start=10_000,
        end=10_400,
        outcome=ms.DECREASE,
        gene="G1",
        cell="K562",
    )
    monkeypatch.setattr(s4.ho, "all_units", lambda: {"crispri:A": (unit,)})
    monkeypatch.setattr(s4.ho, "crispri_references", lambda: {})
    saved = {}
    monkeypatch.setattr(
        s4, "save_result", lambda name, payload, where, manifest: saved.update(n=name, p=payload, m=manifest)
    )
    monkeypatch.setattr(s4, "manifest", lambda programs, params, also=(): {"parameters": params})
    s4.main(["--compiled", str(tmp_path), "--results", str(tmp_path)])
    assert saved["n"] == "attribution_correctness" and "rule" not in saved["p"]
    assert "judge_rule" not in saved["m"]["parameters"]  # a result that names no rule was judged under v1
    (v,) = [v for v in saved["p"]["judged"] if v["axis"] == co.ACTIVITY]
    assert (v["verdict"], v["refutable"]) == (co.INCORRECT, False) and "cross_cell" not in v
