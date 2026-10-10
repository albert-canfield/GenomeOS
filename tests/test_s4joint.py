# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2d: signaltrace's bridge stage re-run with combine2's joint fit (genomeos/attribution/s4joint.py).

Every program here is synthetic, in the compiled generator's own line format, so each joint stage is
checked on the committed code it calls (signaltrace.chromosome_lookup, compile.gate_program,
combine2.gene_rules / classify / fit / run_gene_form, runtime.grn) without the generated store; the
registered run reads that store.
"""

from __future__ import annotations

import copy
import json
import math
import random
from pathlib import Path

import pytest

from genomeos.attribution import combine2 as c2
from genomeos.attribution import crispri
from genomeos.attribution import regulators as regs
from genomeos.attribution import s4joint as sj
from genomeos.attribution import signaltrace as st
from genomeos.genome.reader import PeakIndex
from genomeos.results import RESULTS_DIR


def _element(name: str, start: int, gene: str, effect: str) -> str:
    return (
        f"element {name} {{\n  class: enhancer\n  locus: chrT:{start}-{start + 100}\n  targets: {gene}\n"
        f'  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect'
        f" {effect} log2 fold change on deletion, probability unavailable\n}}\n"
    )


def _rule(name: str, action: str, gene: str, effect: str, cell: str = "K562") -> str:
    strength = min(1.0, abs(float(effect)))
    return (
        f"rule {name} {action} {gene} {{ strength: {strength}; when: cell_type = {cell}; evidence: predicted"
        f' "AlphaGenome deletion, {cell}" effect {effect} log2 fold change on deletion, probability'
        " unavailable }\n"
    )


GENES = ("TGT", "TWO", "BIG", "ONE", "SUP_1")
PROGRAM = (
    "module human.noncoding.chrT\n\n"
    + "".join(
        f'gene {g} {{ symbol: {g.replace("_", "-")}; evidence: curated "GENCODE v50" symbol only, a stub }}\n'
        for g in GENES
    )
    + _element("E1", 100, "TGT", "-0.5")
    + _rule("E1", "activates", "TGT", "-0.5", "placenta")
    + _element("E5", 2000, "TWO", "-0.4")
    + _rule("E5", "activates", "TWO", "-0.4")
    + _element("E6", 2200, "TWO", "-0.2")
    + _rule("E6", "activates", "TWO", "-0.2")
    + _element("E7", 3000, "BIG", "-2.0")
    + _rule("E7", "activates", "BIG", "-2.0")
    + _element("E8", 3200, "BIG", "-2.0")
    + _rule("E8", "activates", "BIG", "-2.0")
    + _element("E9", 4000, "ONE", "-1.0")
    + _rule("E9", "activates", "ONE", "-1.0")
    + _element("E2", 5000, "SUP_1", "+0.3")
    + _rule("E2", "inhibits", "SUP_1", "+0.3")
    + "element E1_measured {\n  class: enhancer\n  locus: chrT:100-200\n  targets: TGT\n"
    '  evidence: experimental "CRISPRi" links TGT in K562 from 1 pairs\n  confidence: 0.90\n}\n'
    'rule E1_measured activates TGT { strength: 0.2; when: cell_type = K562; evidence: experimental "CRISPRi,'
    ' silenced in K562, held-out split, evaluation only, never a feature"; confidence: 0.90 }\n'
)


class FakeReaders:
    """E1, E5, E7 and E9 open in K562; E2, E6 and E8 open nowhere; no peak file for testis."""

    def _peaks(self, biosample: str) -> list[tuple[int, int, float]] | None:
        if biosample == "testis":
            return None
        peaks = [(10, 20, 1.0), (9000, 9100, 1.0)]
        if biosample == "K562":
            peaks += [(150, 160, 5.0), (2050, 2060, 5.0), (3050, 3060, 5.0), (4050, 4060, 5.0)]
        if biosample == "HepG2":
            peaks.append((2250, 2260, 5.0))
        return peaks

    def span(self, biosample: str, chrom: str):
        p = self._peaks(biosample)
        return None if p is None else (min(s for s, _, _ in p), max(e for _, e, _ in p))

    def index(self, biosample: str, chrom: str):
        p = self._peaks(biosample)
        return None if p is None else PeakIndex(p)


CLUSTERS = ["chrT\t2050\t2060\tGATA1\t1000\tK562,HepG2"]  # occupies E5 in K562


def pair(start: int, end: int, gene: str) -> crispri.Pair:
    p = crispri.Pair("chrT", start, end, gene, "K562", "test", 5000.0, 1.0, 1.0, False)
    p.features = {"top_target": 1.0, "deletion_drop": 0.3}
    return p


PAIRS = {
    "two": pair(2050, 2060, "TWO"),  # E5: TWO's open element
    "two6": pair(2250, 2260, "TWO"),  # E6: TWO's closed element
    "big": pair(3050, 3060, "BIG"),  # summed drops 1.5: no positive basal rate
    "one": pair(4050, 4060, "ONE"),  # one rule: the joint fit at k = 1
    "sup": pair(5050, 5060, "SUP-1"),  # inhibitor only
    "tgt": pair(150, 170, "TGT"),  # compiled in placenta, offered in K562; its CRISPRi twin removed
    "none": pair(7000, 7100, "TGT"),  # no compiled element
}
RHO = {"E5": 2**-0.4, "E6": 2**-0.2}


@pytest.fixture(scope="module")
def look() -> st.ChromLookup:
    return st.chromosome_lookup(
        "chrT", PROGRAM, PAIRS.values(), {}, 10.0, None, None,
        regs.parse_clusters(CLUSTERS, ["chrT"]), FakeReaders(), RESULTS_DIR,
    )  # fmt: skip


@pytest.fixture(scope="module")
def joint(look) -> sj.JointLookup:
    return sj.joint_lookup("chrT", PROGRAM, look, {}, 10.0, None, None, FakeReaders(), RESULTS_DIR)


@pytest.fixture(scope="module")
def traced(look, joint) -> dict[str, tuple[dict[str, float], dict[str, str]]]:
    out = {}
    for k, p in PAIRS.items():
        f, why = st.trace_pair(p, look)
        g, why2 = sj.trace_joint(p, look, joint)
        out[k] = ({**f, **g}, {**why, **why2})
    return out


# ---- the registered words --------------------------------------------------------------------------------


def test_the_forms_and_the_decision_are_the_committed_ones_word_for_word():
    assert sj.FORMS == c2.FORMS == ("mean", "competition")
    assert sj.THETA == c2.THETA == 0.5 and sj.FIT is c2.FIT and sj.FORM_SOURCES is c2.FORM_SOURCES
    for mine, theirs in (
        (sj.METRIC, st.METRIC), (sj.BOOTSTRAP, st.BOOTSTRAP), (sj.DECISION, st.DECISION),
        (sj.OUTCOMES, st.OUTCOMES), (sj.CHAIN_OUTCOMES, st.CHAIN_OUTCOMES), (sj.ANCHOR, st.ANCHOR),
        (sj.PREDICTED_ONLY, st.PREDICTED_ONLY), (sj.LOSS, st.LOSS),
    ):  # fmt: skip
        assert mine is theirs
    assert set(sj.LIMITATIONS) >= set(st.LIMITATIONS)
    assert sj.BRIDGE_OF.keys() == set(st.OUTCOMES) - {"HARNESS_DEFECT", "ANCHOR"}
    assert {sj.BRIDGE_OF[k] for k in st.LOSS} == {"LOST"}


def test_each_form_has_its_own_chain_and_every_stage_follows_one_that_exists():
    for form in sj.FORMS:
        chain = sj.CHAINS[form]
        assert chain[:3] == ("s1_anchor", "s2_target", "s3_compiled_value")
        assert chain[3:] == (f"s4_joint_{form}", f"s5_reader_{form}", f"s7_run_{form}")
        assert sj.BRANCHES[form] == (f"s6_regulator_{form}",)
        assert sj.FORM_WORDS[form] in sj.STAGES[f"s4_joint_{form}"]["name"]
    for s in sj.JOINT_STAGES:
        assert sj.STAGES[s]["follows"] in sj.STAGES
    assert json.loads(json.dumps(sj.registered_fields())) == sj.registered_fields()


# ---- the joint fit on the synthetic program --------------------------------------------------------------


def test_every_reached_gene_is_read_and_fitted_as_combine2_fits_it(joint):
    assert joint.reached_genes == {"TGT", "TWO", "BIG", "ONE", "SUP_1"}
    assert set(joint.genes) == joint.reached_genes
    two = joint.genes["TWO"]
    assert two.elements == ["E5", "E6"] and two.k562 == {"E5": "open", "E6": "closed"}
    assert two.rho == pytest.approx(RHO)
    for form in sj.FORMS:
        assert two.forms[form]["ok"] and two.forms[form]["harness_miss"] == []
        assert two.forms[form]["q"] == pytest.approx(c2.fit(form, RHO).q)
        big = joint.genes["BIG"].forms[form]
        assert (big["ok"], big["reason"]) == (False, "basal_not_positive")
        assert big["q"] == pytest.approx(c2.fit(form, {"E7": 0.25, "E8": 0.25}).q)
        assert joint.genes["SUP_1"].forms[form] == {"ok": False, "reason": "no_activating_rule"}
    assert joint.genes["TGT"].elements == ["E1"]  # the experimental twin was removed before gating


def test_the_k562_hold_is_the_forms_closed_form_with_the_other_elements_as_k562_states_them(joint):
    for form in sj.FORMS:
        f = c2.fit(form, RHO)
        want = math.log2(c2.level(f, {"E5": True, "E6": False}) / c2.level(f, {"E5": False, "E6": False}))
        got = joint.genes["TWO"].forms[form]["change"]
        assert set(got) == {"E5"}  # E6 is closed in K562: no hold run
        assert got["E5"] == pytest.approx(want, abs=sj.CLOSED_FORM_TOLERANCE)
        assert got["E5"] != pytest.approx(0.4, abs=1e-3)  # not the observation: E6 is held in K562
        one = joint.genes["ONE"].forms[form]["change"]["E9"]
        assert one == pytest.approx(1.0, abs=sj.CLOSED_FORM_TOLERANCE)  # one element: the observation


# ---- the pairs through the joint stages ------------------------------------------------------------------


def test_a_pair_is_kept_at_the_joint_bridge_where_the_bridge_kept_none(traced):
    f, why = traced["two"]
    assert f["s4_bridge_kept"] == 0.0 and why["s4_bridge"] == "not_identifiable"
    for form in sj.FORMS:
        assert f[f"s4_joint_{form}_kept"] == 1.0
        assert f[f"s4_joint_{form}_drop"] == pytest.approx(0.4)
        assert f[f"s5_reader_{form}_kept"] == 1.0 and f[f"s6_regulator_{form}_kept"] == 1.0
        assert f[f"s7_run_{form}_kept"] == 1.0
        # E6 is held in K562, so removing E5 there is not its own observation (0.4)
        assert f[f"s7_run_{form}_drop"] > 0.0 and f[f"s7_run_{form}_drop"] != pytest.approx(0.4, abs=1e-3)


def test_each_drop_names_its_reason(traced):
    for form in sj.FORMS:
        f, why = traced["two6"]
        assert f[f"s4_joint_{form}_kept"] == 1.0 and f[f"s5_reader_{form}_kept"] == 0.0
        assert why[f"s5_reader_{form}"] == "reader_closed" and why[f"s7_run_{form}"] == "held in K562"
        f, why = traced["big"]
        assert f[f"s4_joint_{form}_kept"] == 0.0 and why[f"s4_joint_{form}"] == "basal_not_positive"
        f, why = traced["sup"]
        assert f[f"s4_joint_{form}_kept"] == 0.0 and why[f"s4_joint_{form}"] == "no_activating_rule"
        f, why = traced["one"]
        assert f[f"s7_run_{form}_drop"] == pytest.approx(1.0, abs=sj.CLOSED_FORM_TOLERANCE)
        assert why[f"s6_regulator_{form}"] == "no factor cluster names K562 at the element"
        f, why = traced["none"]
        assert not any(f[f"{s}_kept"] for s in sj.JOINT_STAGES) and f"s4_joint_{form}" not in why


def test_the_joint_columns_never_read_the_label(look, joint):
    p, q = copy.deepcopy(PAIRS["two"]), copy.deepcopy(PAIRS["two"])
    q.regulated = not p.regulated
    assert sj.trace_joint(p, look, joint) == sj.trace_joint(q, look, joint)


# ---- the harness ---------------------------------------------------------------------------------------


def test_the_harness_passes_a_sound_trace_and_names_each_defect(look, joint, traced):
    rows = []
    for k, p in PAIRS.items():
        p.features.update(traced[k][0])
        rows.append(p)
    assert st.harness({"chrT": look}, rows) == []
    assert sj.harness({"chrT": look}, {"chrT": joint}, rows) == []
    bad = copy.deepcopy(joint)
    bad.genes["TWO"].forms["mean"]["harness_miss"] = ["K562 hold E5: run 1, closed form 2"]
    bad.genes["TWO"].k562["E5"] = "closed"
    del bad.genes["ONE"]
    problems = sj.harness({"chrT": look}, {"chrT": bad}, [])
    assert any("not read from the gated program" in x for x in problems)
    assert any("joint K562 state closed against open" in x for x in problems)
    assert any("K562 hold E5" in x for x in problems)
    stray = pair(1, 2, "X")
    stray.features.update({f"{s}_{c}": 0.0 for s in sj.STAGES for c in ("kept", "drop")})
    stray.features["s7_run_mean_kept"] = 1.0
    assert sj.harness({}, {}, [stray]) == ["s7_run_mean keeps 1 pairs that s5_reader_mean dropped"]


# ---- the readings ----------------------------------------------------------------------------------------


def _row(outcome: str) -> dict:
    return {"outcome": outcome, "reading": f"<{outcome}>"}


def test_the_chain_and_the_bridge_answer_follow_signaltraces_walk():
    stages = {
        "s2_target": _row("WEAKENED"),
        "s3_compiled_value": _row("WEAKENED"),
        "s4_joint_mean": _row("CARRIES"),
        "s5_reader_mean": _row("NONE_KEPT"),
        "s7_run_mean": _row("NONE_KEPT"),
        "s4_joint_competition": _row("UNDECIDED"),
        "s5_reader_competition": _row("NONE_KEPT"),
        "s7_run_competition": _row("NONE_KEPT"),
    }
    kind, text, at = sj.chain_verdict(stages, sj.CHAINS["mean"])
    assert (kind, at) == ("LOST_AT", "s5_reader_mean")
    assert text.endswith("<NONE_KEPT>; it was first weakened at s2_target (target choice)")
    assert sj.bridge_verdict(stages, "mean")[0] == "SURVIVES"
    kind, _text, at = sj.chain_verdict(stages, sj.CHAINS["competition"])
    assert (kind, at) == ("UNDECIDED_AT", "s4_joint_competition")
    b, btext = sj.bridge_verdict(stages, "competition")
    assert b == "UNDECIDED" and "is not decided" in btext and sj.FORM_WORDS["competition"] in btext
    stages["s4_joint_mean"] = _row("WORSE")
    assert sj.bridge_verdict(stages, "mean")[0] == "LOST"


def _population(n_chroms: int = 12, per: int = 12, seed: int = 3) -> list[crispri.Pair]:
    rng = random.Random(seed)
    out = []
    for c in range(n_chroms):
        for i in range(per):
            y = rng.random() < 0.25
            p = crispri.Pair(f"chr{c + 1}", i * 10, i * 10 + 5, f"G{c}_{i}", "K562", "t", 1e4, 1.0, 1.0, y)
            signal = (1.0 if y else 0.0) + rng.gauss(0, 0.6)
            p.features = {
                "log_distance": math.log(1000 + rng.random() * 1e5),
                "log_activity": rng.random(),
                "activity_over_distance": rng.random() - 5,
                "top_target": 1.0,
                "deletion_drop": max(0.0, signal),
            }
            p.covered = True
            out.append(p)
    return out


def test_score_reads_each_form_against_the_anchor_and_the_comparators_reproduce():
    fit, pop = _population(seed=1), _population(seed=2)
    for p in fit + pop:
        f = p.features
        for s in sj.STAGES:
            f[f"{s}_kept"], f[f"{s}_drop"] = 0.0, 0.0
        for s in ("s1_anchor", "s2_target", "s3_compiled_value", "s4_joint_mean"):
            f[f"{s}_kept"], f[f"{s}_drop"] = f["top_target"], f["deletion_drop"]
    anchor = st.anchor_reading(fit, pop)
    rows, scores, weights = sj.comparators(fit, pop, anchor)
    assert sj.comparator_mismatch(rows, copy.deepcopy(rows)) == []
    moved = copy.deepcopy(rows)
    moved["s3_compiled_value"]["gain_over_activity_distance"]["gain"] += 0.0001
    problems = sj.comparator_mismatch(rows, moved)
    assert len(problems) == 1 and problems[0].startswith("s3_compiled_value gain_over_activity_distance gain")
    out = sj.score_forms(fit, pop, anchor, rows, scores, weights)
    s = out["stages"]
    assert s["s4_joint_mean"]["gain_over_activity_distance"] == s["s2_target"]["gain_over_activity_distance"]
    assert s["s4_joint_mean"]["against_the_stage_it_follows"]["gain"] == 0.0
    assert s["s4_joint_competition"]["outcome"] == "NONE_KEPT"
    comp = out["per_form"]["competition"]
    assert (comp["bridge"], comp["chain"]) == ("LOST", "LOST_AT")
    assert comp["stage_where_it_stops"] in ("s2_target", "s3_compiled_value", "s4_joint_competition")
    assert out["per_form"]["mean"]["bridge"] == sj.BRIDGE_OF[s["s4_joint_mean"]["outcome"]]
    assert out["reading"].startswith("with the genes fitted jointly under the runtime's mean")
    assert comp["reading"] in out["reading"]


# ---- the committed inputs this lane binds to ------------------------------------------------------------


def test_combine2s_committed_registration_is_combine2_py_on_every_field_the_run_reads():
    import tracked_paths as tp

    path = Path("data/results/combine2_registration.json")
    tp.must_be_present(path)
    assert sj.combine2_mismatch(json.loads(path.read_text())) == []
    assert sj.combine2_mismatch({**json.loads(path.read_text()), "theta": 0.4}) == ["theta"]
