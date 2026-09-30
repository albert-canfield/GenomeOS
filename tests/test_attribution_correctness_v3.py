# SPDX-License-Identifier: AGPL-3.0-or-later
"""The target rule, versioned (registered 2026-09-29, lane-judge3): v3 judges a target claim only in the
cell it states, for supported and refuted verdicts alike, and keeps every other cell beside the verdict;
v1 and v2 stay selectable and still reproduce their committed files."""

from __future__ import annotations

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
V3_RESULT = ROOT / "data/results/attribution_correctness_v3.json"
E = ("chr1", 10_000, 10_400)
CELLS = ("K562", "HepG2", "placenta", "", "unknown")
KINDS = ("crispri_decrease", "crispri_increase", "crispri_null_well_powered", "crispri_null_underpowered")
RESPONSES = ("crispri_decrease", "crispri_increase")
ELSEWHERE_READ = (*RESPONSES, "crispri_null_well_powered")  # what a target's cross-cell finding keeps
V1_NULL_TO_CONTEXT = "responds in another cell; the stated cell's null is read on context"


def claim(axis: str, value: str, gene: str = "G1", cell: str = "K562") -> co.Claim:
    return co.Claim("EH1", *E, axis, value, gene, cell)


def obs(kind: str, cell: str = "K562", gene: str = "G1") -> co.Observation:
    return co.Observation(kind, f"crispri:{cell or 'x'}", gene, cell, "training")


def target(cell: str = "K562") -> co.Claim:
    return claim(co.TARGET, "G1", cell=cell)


def in_cell(seen: list[co.Observation], cell: str) -> list[co.Observation]:
    return [o for o in seen if co.stated(cell) and co.norm_cell(o.cell) == co.norm_cell(cell)]


# --- the rule ----------------------------------------------------------------------------------------
def test_v3_is_the_default_and_v1_and_v2_stay_selectable():
    assert co.DEFAULT_RULE == co.RULE_V3 and co.RULES == (co.RULE_V1, co.RULE_V2, co.RULE_V3)
    seen = [obs("crispri_decrease", "K562")]
    assert co.verdict_of(target("placenta"), seen).rule == co.RULE_V3
    assert [co.verdict_of(target("placenta"), seen, r).verdict for r in co.RULES] == [
        co.CORRECT,
        co.CORRECT,
        co.NOT_JUDGED,
    ]
    with pytest.raises(ValueError):
        co.verdict_of(target(), [], rule="v4")
    with pytest.raises(ValueError):
        co.judge([], ho.Labels("x", lambda u, e: None, frozenset()), units={}, sources=[], rule="v4")


@pytest.mark.parametrize("cell", CELLS)
def test_v3_never_decides_a_target_from_another_cell(cell):
    """Every combination of up to three CRISPRi outcomes in two cells: v3's target verdict is what the
    stated cell's observations alone give, so another cell neither establishes nor refutes the target, and
    neither enables nor blocks a refutation; each response and well-powered null elsewhere is kept beside
    the verdict, and a claim the stated cell does not decide is not assessed whenever another cell
    responded."""
    outcomes = [(k, c) for k in KINDS for c in ("K562", "HepG2")]
    c = target(cell)
    for n in range(4):
        for combo in itertools.combinations(outcomes, n):
            seen = [obs(k, oc) for k, oc in combo]
            here = in_cell(seen, cell)
            v, alone, v1 = co.verdict_of(c, seen), co.verdict_of(c, here), co.verdict_of(c, seen, co.RULE_V1)
            assert v.rule == co.RULE_V3 and all(o in here for o in v.deciding), combo
            if v.verdict in co.JUDGED or alone.verdict in co.JUDGED:
                assert (v.verdict, v.detail, v.deciding, v.refutable) == (
                    alone.verdict,
                    alone.detail,
                    alone.deciding,
                    True,
                ), combo
            elif any(o.kind in RESPONSES for o in seen):  # the stated cell decides nothing, another responded
                assert (v.reason, v.refutable, v.deciding) == (co.NOT_ASSESSED, False, ()), combo
            else:  # no cell decides it under any rule: v1's reason and detail
                assert (v.verdict, v.reason, v.detail, v.cross_cell) == (v1.verdict, v1.reason, v1.detail, ())
            if v.verdict in co.JUDGED or v.reason == co.NOT_ASSESSED:
                kept = [o for o in seen if o not in here and o.kind in ELSEWHERE_READ]
                assert [o for o, _ in v.cross_cell] == kept, combo
                for o, f in v.cross_cell:
                    assert f == (co.AGREES if o.kind in RESPONSES else co.DISAGREES)
            if v1.verdict in co.JUDGED and alone.verdict not in co.JUDGED:  # v1 decided it only elsewhere
                assert v.reason == co.NOT_ASSESSED, combo


def test_the_v1_condition_on_a_refutation_is_not_kept_and_what_follows():
    """TARGET_V1_CONDITION: a response elsewhere no longer turns a stated-cell null into a context error;
    the null refutes the target in context, and the unchanged context rule refutes the context claim."""
    seen = [obs("crispri_decrease", "HepG2"), obs("crispri_null_well_powered", "K562")]
    t1, t3 = co.verdict_of(target(), seen, co.RULE_V1), co.verdict_of(target(), seen)
    assert (t1.verdict, t1.detail) == (co.CORRECT, V1_NULL_TO_CONTEXT)
    assert (t3.verdict, t3.detail, t3.deciding) == (
        co.INCORRECT,
        "well-powered null in the stated cell",
        (seen[1],),
    )
    assert t3.cross_cell == ((seen[0], co.AGREES),)  # the response elsewhere, kept as supporting evidence
    ctx = claim(co.CONTEXT, "K562")
    assert co.verdict_of(ctx, seen, co.RULE_V1).verdict == co.verdict_of(ctx, seen).verdict == co.INCORRECT
    # a null in the stated cell refutes by itself: no response in any cell is needed
    alone = [obs("crispri_null_well_powered", "K562")]
    assert co.verdict_of(target(), alone).verdict == co.INCORRECT
    assert co.verdict_of(ctx, alone).reason == co.CONDITION_UNMET  # context's own condition, unchanged
    # the condition S4 does put on a null, its power, stays: an underpowered null decides nothing
    weak = [obs("crispri_null_underpowered", "K562"), obs("crispri_decrease", "HepG2")]
    assert co.verdict_of(target(), weak).reason == co.NOT_ASSESSED
    assert co.verdict_of(target(), weak[:1]).reason == co.CANNOT_ONLY


def test_a_claim_with_no_stated_cell_is_never_judged_on_target():
    for cell in ("", "unknown"):
        v = co.verdict_of(target(cell), [obs("crispri_decrease")])
        assert (v.verdict, v.reason, v.detail) == (co.NOT_JUDGED, co.NOT_ASSESSED, "the claim states no cell")
        v = co.verdict_of(target(cell), [obs("crispri_null_well_powered")])
        assert (v.verdict, v.reason, v.detail) == (co.NOT_JUDGED, co.OUTSIDE_SCOPE, co.NULL_ELSEWHERE)


def test_only_the_target_moves_between_v2_and_v3():
    seen = [
        obs("crispri_decrease", "K562"),
        obs("crispri_increase", "HepG2"),
        obs("crispri_null_well_powered", "WTC11"),
        co.Observation("lentimpra_active", "lentimpra:K562", cell="K562"),
        co.Observation("vista_positive", "vista", cell=co.VISTA_CONTEXT),
        co.Observation("registry_biochemical", "encode", value="enhancer_like"),
        co.Observation("sequence_annotation", "rmsk", value="unique"),
    ]
    others = [
        claim(co.ACTIVITY, value, cell=cell)
        for value in (*co.DIRECTION, "no_effect_measured")
        for cell in ("K562", "HepG2", "WTC11", "placenta", "")
    ]
    others += [claim(co.CONTEXT, cell, cell=cell) for cell in ("K562", "HepG2", "WTC11", "placenta")]
    others += [
        claim(co.ACTIVITY, "active_in_reporter", gene="", cell="HepG2"),
        claim(co.ACTIVITY, "inactive_in_reporter", gene="", cell=""),
        claim(co.ROLE, "enhancer_like", gene="", cell=""),
        claim(co.ORIGIN, "unique", gene="", cell=""),
    ]
    for c in others:
        v2, v3 = co.verdict_of(c, seen, co.RULE_V2), co.verdict_of(c, seen)
        assert s4._pair(v2) == s4._pair(v3), c
        assert (v2.refutable, v2.cross_cell) == (v3.refutable, v3.cross_cell), c
    # on target, v2 decided every one of these with another cell's response; v3 with the stated cell's only
    got = {cell: co.verdict_of(target(cell), seen) for cell in ("K562", "HepG2", "WTC11", "placenta", "")}
    assert {cell: (v.verdict, v.reason) for cell, v in got.items()} == {
        "K562": (co.CORRECT, ""),
        "HepG2": (co.CORRECT, ""),
        "WTC11": (co.INCORRECT, ""),
        "placenta": (co.NOT_JUDGED, co.NOT_ASSESSED),
        "": (co.NOT_JUDGED, co.NOT_ASSESSED),
    }
    assert (got["K562"].deciding, got["HepG2"].deciding) == ((seen[0],), (seen[1],))
    assert all(co.verdict_of(target(cell), seen, co.RULE_V2).verdict == co.CORRECT for cell in got)


def test_the_elsewhere_evidence_survives_serialisation():
    seen = [obs("crispri_decrease", "K562"), obs("crispri_null_well_powered", "HepG2")]
    v = co.verdict_of(target("Whole_Blood"), seen)
    assert (v.reason, [f for _, f in v.cross_cell]) == (co.NOT_ASSESSED, [co.AGREES, co.DISAGREES])
    d = json.loads(json.dumps(v.to_dict()))
    assert d["rule"] == co.RULE_V3 and [x["finding"] for x in d["cross_cell"]] == [co.AGREES, co.DISAGREES]
    assert co.Verdict.from_dict(d) == v
    kept_v1 = co.verdict_of(target("Whole_Blood"), seen[1:])
    assert (
        kept_v1.cross_cell == ()
        and co.Verdict.from_dict(json.loads(json.dumps(kept_v1.to_dict()))) == kept_v1
    )


def _units() -> dict[str, tuple[ho.Unit, ...]]:
    def unit(start: int, outcome: str, gene: str, cell: str) -> ho.Unit:
        return ho.Unit(
            source="crispri:A",
            chrom="chr1",
            start=start,
            end=start + 400,
            outcome=outcome,
            gene=gene,
            cell=cell,
        )

    return {
        "crispri:A": (
            unit(10_000, ms.DECREASE, "G1", "K562"),  # EH1: G1 fell in K562; the claim states Whole_Blood
            unit(20_000, ms.DECREASE, "G2", "K562"),  # EH2: G2 fell in K562, the cell the claim states
            unit(20_000, ms.NULL_INFORMATIVE, "G2", "HepG2"),  # ... and a null elsewhere, kept beside
            unit(30_000, ms.NULL_INFORMATIVE, "G3", "K562"),  # EH3: a null in the stated cell only
        )
    }


def _claims() -> list[co.Claim]:
    def rule(n: int, start: int, gene: str, cell: str, act: str) -> list[co.Claim]:
        c = ("chr1", start, start + 400)
        return [
            co.Claim(f"EH{n}", *c, co.TARGET, gene, gene, cell),
            co.Claim(f"EH{n}", *c, co.ACTIVITY, act, gene, cell),
            co.Claim(f"EH{n}", *c, co.CONTEXT, cell, gene, cell),
        ]

    return [
        *rule(1, 10_000, "G1", "Whole_Blood", "represses_target"),
        *rule(2, 20_000, "G2", "K562", "activates_target"),
        *rule(3, 30_000, "G3", "K562", "activates_target"),
    ]


def _report(rule: str = co.RULE_V3) -> co.Report:
    labels = ho.Labels("synthetic", lambda u, e: None, frozenset({ho.MODEL}))
    return co.judge(_claims(), labels, units=_units(), sources=["crispri:A"], rule=rule)


def test_the_report_counts_in_the_owners_terms_and_round_trips():
    r = _report()
    ic = r.to_dict()["in_context"]
    assert ic[co.TARGET] == {
        "supported_in_context": 1,
        "refuted_in_context": 1,
        "unresolved_in_context": 0,
        "observation_model_inadequate_in_context": 0,
        "unassessed_in_context": 1,
        "supported_elsewhere_only": 1,
        "supported_somewhere": 2,
    }
    assert (ic[co.ACTIVITY]["supported_in_context"], ic[co.ACTIVITY]["unassessed_in_context"]) == (1, 1)
    assert ic[co.ACTIVITY]["supported_somewhere"] == 1  # EH1's K562 fall disagrees with its repression
    assert r.axes[co.TARGET].to_dict()["cross_cell_findings"] == {
        co.CORRECT: {co.DISAGREES: 1},
        co.NOT_ASSESSED: {co.AGREES: 1},
    }
    d = json.loads(json.dumps(r.to_dict()))
    back = co.Report.from_dict(d)
    assert back.to_dict() == d and back.rule == co.RULE_V3
    assert [(v.claim.axis, v.claim.element) for v in back.not_assessed] == [
        (co.TARGET, "EH1"),
        (co.ACTIVITY, "EH1"),
    ]


def _keys_and_strings(x, keys: set, strings: list) -> None:
    if isinstance(x, dict):
        for k, v in x.items():
            keys.add(str(k))
            _keys_and_strings(v, keys, strings)
    elif isinstance(x, list):
        for v in x:
            _keys_and_strings(v, keys, strings)
    elif isinstance(x, str):
        strings.append(x)


def assert_never_called_accuracy(block) -> None:
    keys, strings = set(), []
    _keys_and_strings(block, keys, strings)
    assert not [k for k in keys if "accuracy" in k.lower()]
    assert all(s == co.IN_CONTEXT_CAUTION for s in strings if "accuracy" in s.lower())
    assert co.IN_CONTEXT_CAUTION in strings


def test_the_report_never_labels_these_counts_accuracy():
    assert "not general accuracy" in co.IN_CONTEXT_CAUTION and "39 of 39" in co.IN_CONTEXT_CAUTION
    for text in (co.IN_CONTEXT_SCOPE, *co.IN_CONTEXT_TERMS, *co.IN_CONTEXT_TERMS.values()):
        assert "accuracy" not in text.lower()
    d = _report().to_dict()
    assert_never_called_accuracy(d["in_context"])
    for block, quantity, axis in co.CAUTIONED:  # every established-over-decided share S4 prints
        assert d[block][quantity][axis]["caution"] == co.IN_CONTEXT_CAUTION
    for rule in (co.RULE_V1, co.RULE_V2):  # the earlier reports are written as they were
        e = _report(rule).to_dict()
        assert "in_context" not in e
        assert all("caution" not in e[b][q][a] for b, q, a in co.CAUTIONED)
    assert_never_called_accuracy(s4.in_context_side_by_side({r: _report(r) for r in co.RULES}))


def test_v2s_rule_record_and_s4s_registration_are_what_their_files_hold():
    assert co.rule_registration()["rules"] == [co.RULE_V1, co.RULE_V2]
    assert co.rule_registration()["default"] == co.RULE_V2
    v2 = json.loads(V2_RESULT.read_text())
    assert json.loads(json.dumps(co.rule_registration())) == v2["rule_registration"]
    assert v2["result_manifest"]["parameters"]["rule_registration"] == v2["rule_registration"]
    v1 = json.loads(V1_RESULT.read_text())
    assert json.loads(json.dumps(co.registration())) == v1["registration"]
    r3 = co.rule_registration_v3()
    assert (r3["rules"], r3["default"], r3["target_rule"][co.RULE_V3]) == (
        list(co.RULES),
        co.RULE_V3,
        co.TARGET_RULE[co.RULE_V3],
    )


def test_the_s4_script_keeps_its_names_and_writes_v3_beside():
    assert s4.NAMES == {co.RULE_V1: "attribution_correctness", co.RULE_V2: "attribution_correctness_v2"}
    assert s4.NAME_V3 == "attribution_correctness_v3"
    src = (ROOT / "scripts/s4_correctness_run.py").read_text()
    assert "report = judged_under(co.RULE_V1, a, units, refs)" in src  # no --rule: v1, as S4 committed


def test_the_document_carries_the_registration_word_for_word():
    doc = (ROOT / "docs/ATTRIBUTION.md").read_text()
    start = doc.index("## The target rule, versioned")
    section = re.sub(r"\s+", " ", re.sub(r"\n> ?", "\n", doc[start : doc.index("\n## ", start + 5)]))
    for text in (
        co.RULE_DECISION_V3.split(": ", 1)[1],
        co.TARGET_RULE[co.RULE_V3],
        co.TARGET_V1_CONDITION,
        co.IN_CONTEXT_SCOPE,
        *co.IN_CONTEXT_TERMS.values(),
        co.IN_CONTEXT_CAUTION,
        co.RULE_RECORD_V3,
        co.SEEN_BEFORE_V3,
        co.EXPECTED_V3,
    ):
        assert re.sub(r"\s+", " ", text) in section
    assert co.TARGET_RULE[co.RULE_V1].endswith(co.AXIS_RULE[co.TARGET])


# --- the v3 run, end to end on a synthetic program ---------------------------------------------------
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


def test_the_v3_run_reports_three_rules_and_every_moved_claim(tmp_path, monkeypatch):
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
        s4,
        "save_result",
        lambda name, payload, where, manifest: saved.update(name=name, p=payload, m=manifest),
    )
    monkeypatch.setattr(s4, "manifest", lambda programs, params, also=(): {"parameters": params})
    s4.main(["--rule", "v3", "--compiled", str(tmp_path), "--results", str(tmp_path)])
    p = saved["p"]
    assert saved["name"] == "attribution_correctness_v3" and p["rule"] == co.RULE_V3
    assert saved["m"]["parameters"]["judge_rule"] == co.RULE_V3
    assert {r: e["rerun_reproduces_it"] for r, e in p["earlier_results"].items()} == {
        co.RULE_V1: None,  # not the registered run: not compared
        co.RULE_V2: None,
    }
    sbs = p["side_by_side"][co.TARGET]
    assert (sbs[co.RULE_V1][co.CORRECT], sbs[co.RULE_V2][co.CORRECT], sbs[co.RULE_V3][co.CORRECT]) == (
        2,
        2,
        1,
    )
    assert sbs[co.RULE_V3]["not_judged_by_reason"][co.NOT_ASSESSED] == 1
    assert p["side_by_side"][co.ACTIVITY][co.RULE_V2] == p["side_by_side"][co.ACTIVITY][co.RULE_V3]
    m = p["moved"]
    assert (m["claims"], m["by_axis"], m["unexplained_count_changes"]) == (1, {co.TARGET: 1}, [])
    (e,) = m["each"]
    assert e["claim"]["element"] == "EH1" and e[co.RULE_V3]["cross_cell"][0]["finding"] == co.AGREES
    assert p["moved_from_v1"]["by_axis"] == {co.TARGET: 1, co.ACTIVITY: 1}
    ic = p["in_context_side_by_side"]
    assert ic[co.RULE_V1][co.TARGET]["decided_with_another_cell"] == 1
    assert (
        ic[co.RULE_V3][co.TARGET]["supported_somewhere"]
        == ic[co.RULE_V1][co.TARGET]["supported_somewhere"]
        == 2
    )
    assert p["in_context"][co.TARGET]["unassessed_in_context"] == 1
    assert_never_called_accuracy(p["in_context"])
    assert_never_called_accuracy(ic)


# --- the committed re-judge --------------------------------------------------------------------------
needs_v3 = pytest.mark.skipif(
    not V3_RESULT.exists(), reason="the v3 re-judge is written once, after the code"
)


@pytest.fixture(scope="module")
def v3() -> dict:
    return json.loads(V3_RESULT.read_text())


@needs_v3
def test_the_v3_result_names_its_rule_and_its_stamp_is_clean(v3):
    m = v3["result_manifest"]
    assert v3["rule"] == co.RULE_V3 and m["parameters"]["judge_rule"] == co.RULE_V3
    assert m["complete"] and m["code"]["dirty"] is False
    assert v3["expected_before_the_run"] == co.EXPECTED_V3 and v3["alphagenome_requests"] == 0


@needs_v3
def test_the_earlier_results_are_kept_unchanged_beside_it(v3):
    read = {i["path"]: i["sha256"] for i in v3["result_manifest"]["inputs"]}
    for rule, path in ((co.RULE_V1, V1_RESULT), (co.RULE_V2, V2_RESULT)):
        e = v3["earlier_results"][rule]
        assert e["rerun_reproduces_it"] is True and e["rule"] == rule
        assert hashlib.sha256(path.read_bytes()).hexdigest() == read[f"data/results/{path.name}"]


@needs_v3
def test_only_the_59_targets_move_between_v2_and_v3(v3):
    sbs = v3["side_by_side"]
    for axis in (co.ORIGIN, co.ROLE, co.ACTIVITY, co.CONTEXT):
        assert sbs[axis][co.RULE_V2] == sbs[axis][co.RULE_V3], axis
    t2, t3 = sbs[co.TARGET][co.RULE_V2], sbs[co.TARGET][co.RULE_V3]
    assert (t3["judged"], t3[co.CORRECT], t3[co.INCORRECT], t3[co.UNRESOLVED]) == (39, 39, 0, 0)
    assert t3["not_judged_by_reason"][co.NOT_ASSESSED] == 59 and t3[co.NOT_JUDGED] - t2[co.NOT_JUDGED] == 59
    for r in co.REASONS:
        assert t3["not_judged_by_reason"][r] == t2["not_judged_by_reason"][r]
    assert (t2["refutable"], t3["refutable"]) == (39, 39)
    m = v3["moved"]
    assert (m["claims"], m["by_axis"], m["unexplained_count_changes"]) == (59, {co.TARGET: 59}, [])
    assert (m["v2_verdicts"], m["v2_decided_in"]) == ({co.CORRECT: 59}, {"another cell only": 59})
    assert m["v3_verdicts"] == {f"{co.NOT_JUDGED}:{co.NOT_ASSESSED}": 59}
    m1 = v3["moved_from_v1"]
    assert (m1["claims"], m1["by_axis"], m1["unexplained_count_changes"]) == (
        118,
        {co.TARGET: 59, co.ACTIVITY: 59},
        [],
    )
    t = v3["axes"][co.TARGET]
    assert t["decided_by_kind"] == {"crispri_decrease": 39} and t["decided_with_heldout_file_pair"] == 1


@needs_v3
def test_the_owners_terms_are_39_59_98_on_target_and_as_v2_on_direction(v3):
    ic = v3["in_context"]
    assert ic[co.TARGET] == {
        "supported_in_context": 39,
        "refuted_in_context": 0,
        "unresolved_in_context": 0,
        "observation_model_inadequate_in_context": 0,
        "unassessed_in_context": 59,
        "supported_elsewhere_only": 59,
        "supported_somewhere": 98,
    }
    assert ic[co.ACTIVITY] == {
        "supported_in_context": 39,
        "refuted_in_context": 0,
        "unresolved_in_context": 0,
        "observation_model_inadequate_in_context": 0,
        "unassessed_in_context": 59,
        "supported_elsewhere_only": 54,
        "supported_somewhere": 93,
    }
    sbs = v3["in_context_side_by_side"]
    assert sbs[co.RULE_V1][co.TARGET]["decided_with_another_cell"] == 59
    assert sbs[co.RULE_V2][co.ACTIVITY] == sbs[co.RULE_V3][co.ACTIVITY]
    assert_never_called_accuracy(ic)
    assert_never_called_accuracy(sbs)


@needs_v3
def test_the_elsewhere_evidence_is_kept_on_every_moved_target(v3):
    for e in v3["moved"]["each"]:
        cross = e[co.RULE_V3]["cross_cell"]
        agrees = {(o["kind"], o["cell"], o["source"]) for o in cross if o["finding"] == co.AGREES}
        assert agrees and {(o["kind"], o["cell"], o["source"]) for o in e[co.RULE_V2]["deciding"]} <= agrees
    assert sum(1 for v in v3["not_assessed"] if v["axis"] == co.TARGET) == 59


# --- v1 and v2 reproduce their committed files, and v3 its re-judge (local data only) ----------------
live = pytest.mark.skipif(
    not (ROOT / ho.COMPILED_DIR).exists() or not (ROOT / "data/knowledge/crispri").exists(),
    reason="the compiled programs and the CRISPRi benchmark are local and untracked",
)


@pytest.fixture(scope="module")
def rerun() -> dict:
    """The v2 run as scripts/s4_correctness_run.py --rule v2 makes it, captured instead of written: it
    reruns v1 and stops unless v1 reproduces attribution_correctness.json; and, once the v3 re-judge is
    committed, v3 judged again on the same claims and sources."""
    out: dict = {}
    with pytest.MonkeyPatch.context() as mp:
        mp.chdir(ROOT)
        mp.setattr(s4, "save_result", lambda name, payload, where, manifest: out.update(name=name, p=payload))
        mp.setattr(s4, "manifest", lambda programs, params, also=(): {"parameters": params})
        s4.main(["--rule", "v2"])
        if V3_RESULT.exists():
            units, refs = ho.all_units(), ho.crispri_references()
            a = Namespace(compiled=ho.COMPILED_DIR, chroms=None)
            out[co.RULE_V3] = s4.judged_under(co.RULE_V3, a, units, refs)
    return out


@live
def test_v1_and_v2_reproduce_their_committed_files_exactly(rerun):
    v2 = json.loads(V2_RESULT.read_text())
    got = json.loads(json.dumps(rerun["p"], default=str))
    assert rerun["name"] == "attribution_correctness_v2"
    assert got["historical_result"]["v1_rerun_reproduces_it"] is True  # v1: the whole file, key by key
    drop = set(s4.RUN_KEYS)
    assert s4._text(got, drop) == s4._text(v2, drop)


@live
@needs_v3
def test_v3_reproduces_the_committed_re_judge(rerun, v3):
    body = json.loads(json.dumps(rerun[co.RULE_V3].to_dict(), default=str))
    for k in (
        "rule",
        "in_context",
        "quantities",
        "also_reported",
        "axes",
        "not_claims",
        "judged",
        "not_assessed",
    ):
        assert body[k] == v3[k], k
