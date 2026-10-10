# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2e: the bridge stage under the multiplicative composition (genomeos/attribution/saturate.py).

Every program here is synthetic, in the compiled generator's own line format, so each stage is checked on
the committed code it calls (signaltrace.chromosome_lookup, compile.gate_program, combine2.gene_rules /
classify, s4joint's records and comparators, runtime.grn's opt-in `multiplicative`) without the generated
store; the registered run reads that store.
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
from genomeos.attribution import saturate as sa
from genomeos.attribution import signaltrace as st
from genomeos.genome.reader import PeakIndex
from genomeos.results import RESULTS_DIR
from genomeos.runtime import grn


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


GENES = ("TGT", "TWO", "BIG", "ONE", "SUP_1", "INH", "MIX")
ELEMENTS = (  # name, start, gene, log2 fold on deletion, action
    ("E5", 2000, "TWO", "-0.4", "activates"),
    ("E6", 2200, "TWO", "-0.2", "activates"),
    ("E7", 3000, "BIG", "-2.0", "activates"),
    ("E8", 3200, "BIG", "-2.0", "activates"),
    ("E9", 4000, "ONE", "-1.0", "activates"),
    ("E2", 5000, "SUP_1", "+0.3", "inhibits"),
    ("E3", 6000, "INH", "-1.0", "activates"),
    ("E4", 6200, "INH", "+1.5", "inhibits"),  # RHO 2.83 > 2: no inhibitor strength in (0, 1]
    ("E10", 8000, "MIX", "-1.5", "activates"),
    ("E11", 8200, "MIX", "+0.5", "inhibits"),
)
PROGRAM = (
    "module human.noncoding.chrT\n\n"
    + "".join(
        f'gene {g} {{ symbol: {g.replace("_", "-")}; evidence: curated "GENCODE v50" symbol only, a stub }}\n'
        for g in GENES
    )
    + _element("E1", 100, "TGT", "-0.5")
    + _rule("E1", "activates", "TGT", "-0.5", "placenta")
    + "".join(_element(n, s, g, e) + _rule(n, a, g, e) for n, s, g, e, a in ELEMENTS)
    + "element E1_measured {\n  class: enhancer\n  locus: chrT:100-200\n  targets: TGT\n"
    '  evidence: experimental "CRISPRi" links TGT in K562 from 1 pairs\n  confidence: 0.90\n}\n'
    'rule E1_measured activates TGT { strength: 0.2; when: cell_type = K562; evidence: experimental "CRISPRi,'
    ' silenced in K562, held-out split, evaluation only, never a feature"; confidence: 0.90 }\n'
)


class FakeReaders:
    """E1, E5, E7, E9, E3, E10 and E11 open in K562; E2, E4, E6 and E8 closed there; no testis file."""

    def _peaks(self, biosample: str) -> list[tuple[int, int, float]] | None:
        if biosample == "testis":
            return None
        peaks = [(10, 20, 1.0), (9000, 9100, 1.0)]
        if biosample == "K562":
            peaks += [(s + 50, s + 60, 5.0) for s in (100, 2000, 3000, 4000, 6000, 8000, 8200)]
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
    "two": pair(2050, 2060, "TWO"),  # E5: TWO's open element, E6 closed in K562
    "two6": pair(2250, 2260, "TWO"),  # E6: TWO's closed element
    "big": pair(3050, 3060, "BIG"),  # summed drops 1.5: the mean and competition cannot fit it
    "one": pair(4050, 4060, "ONE"),  # one rule
    "sup": pair(5050, 5060, "SUP-1"),  # inhibitor only
    "inh": pair(6050, 6060, "INH"),  # an inhibitor with RHO > 2
    "mix": pair(8050, 8060, "MIX"),  # an activator and a feasible inhibitor, both open
    "tgt": pair(150, 170, "TGT"),  # compiled in placenta, offered in K562; its CRISPRi twin removed
    "none": pair(7000, 7100, "TGT"),  # no compiled element
}


@pytest.fixture(scope="module")
def look() -> st.ChromLookup:
    return st.chromosome_lookup(
        "chrT", PROGRAM, PAIRS.values(), {}, 10.0, None, None,
        regs.parse_clusters(CLUSTERS, ["chrT"]), FakeReaders(), RESULTS_DIR,
    )  # fmt: skip


@pytest.fixture(scope="module")
def joint(look) -> sj.JointLookup:
    return sa.joint_lookup("chrT", PROGRAM, look, {}, 10.0, None, None, FakeReaders(), RESULTS_DIR)


@pytest.fixture(scope="module")
def traced(look, joint) -> dict[str, tuple[dict[str, float], dict[str, str]]]:
    out = {}
    for k, p in PAIRS.items():
        f, why = st.trace_pair(p, look)
        g, why2 = sa.trace(p, look, joint)
        out[k] = ({**f, **g}, {**why, **why2})
    return out


# ---- the registered words --------------------------------------------------------------------------------


def test_the_metric_decision_and_sentences_are_m1c_and_m2d_word_for_word():
    for mine, theirs in (
        (sa.METRIC, st.METRIC), (sa.BOOTSTRAP, st.BOOTSTRAP), (sa.DECISION, st.DECISION),
        (sa.OUTCOMES, st.OUTCOMES), (sa.CHAIN_OUTCOMES, st.CHAIN_OUTCOMES), (sa.ANCHOR, st.ANCHOR),
        (sa.PREDICTED_ONLY, st.PREDICTED_ONLY), (sa.LOSS, st.LOSS), (sa.POPULATION, sj.POPULATION),
        (sa.BRIDGE_OF, sj.BRIDGE_OF), (sa.BRIDGE_OUTCOMES, sj.BRIDGE_OUTCOMES),
        (sa.FORM_READING, sj.FORM_READING), (sa.COMPARATORS, sj.COMPARATORS),
    ):  # fmt: skip
        assert mine is theirs
    assert set(sa.LIMITATIONS) >= set(st.LIMITATIONS)
    assert sa.FORM in grn.OPT_IN_ACTIVATOR_RULES and sa.FORM not in grn.ACTIVATOR_RULES
    assert sa.FIT["principle"] is c2.FIT["principle"]
    assert sa.CHAIN == ("s1_anchor", "s2_target", "s3_compiled_value", sa.S4, sa.S5, sa.S7)
    for s in sa.MULT_STAGES:
        assert sa.STAGES[s]["follows"] in sa.STAGES
    assert json.loads(json.dumps(sa.registered_fields())) == sa.registered_fields()


@pytest.mark.parametrize(
    "f",
    [lambda w: w / (1.0 + w), lambda w: w / (0.1 + w), math.sqrt, lambda w: min(w, 0.7), lambda w: w],
)
def test_why_not_saturating_a_concave_composition_never_sums_its_drops_past_the_whole(f):
    """WHY_NOT_SATURATING, numerically: for concave nondecreasing f with f(0) = 0 the single-removal
    drops f(W) - f(W - w_i) add up to at most f(W) - f(0), whatever the elements' sizes."""
    rng = random.Random(7)
    for _ in range(300):
        w = [rng.expovariate(1.0) * rng.choice((0.01, 1.0, 20.0)) for _ in range(rng.randint(2, 30))]
        total = sum(w)
        drops = sum(f(total) - f(total - x) for x in w)
        assert drops <= f(total) - f(0.0) + 1e-12


# ---- the multiplicative fit ------------------------------------------------------------------------------


def test_the_fit_reproduces_every_single_removal_where_the_additive_fit_has_no_basal_rate():
    rho = {"a": 0.25, "b": 0.25, "c": 0.5, "r": 2**0.5}
    assert c2.fit("mean", rho).reason == "basal_not_positive"  # summed drops 2.0
    f = sa.fit(rho)
    assert f.ok and f.summed_drop == pytest.approx(2.0)
    every = dict.fromkeys(rho, True)
    for el, r in rho.items():
        assert sa.level(f, {**every, el: False}) / sa.level(f, every) == pytest.approx(r)
    assert sa.level(f, every) == 1.0
    assert max(f.s[e] for e in f.beta) == 1.0 and all(0 < f.s[e] <= 1 for e in rho)
    assert f.gain * f.s["a"] * sa.H == pytest.approx(-math.log(0.25))
    assert f.s["r"] == pytest.approx((1 - 2**-0.5) / sa.H)
    assert f.vmax_frac == pytest.approx(0.25 * 0.25 * 0.5 / (1 - f.s["r"] * sa.H))


def test_the_fit_refuses_what_the_runtime_cannot_carry():
    assert sa.fit({"a": 0.5, "r": 2.5}).reason == "inhibitor_strength_above_1"
    assert sa.fit({"a": 0.5, "r": 2.0}).ok  # s = 1 exactly
    assert sa.fit({"r": 1.5}).reason == "no_activating_rule"
    assert sa.fit({"a": 1.0, "b": 0.5}).reason == "observation_null_effect"
    deep = {f"e{i}": math.exp(-300.0) for i in range(3)}
    assert sa.fit(deep).reason == "exponent_beyond_float_range"


# ---- the joint lookup on the synthetic program -------------------------------------------------------------


def test_every_reached_gene_is_read_with_combine2s_exclusions_and_fitted(joint):
    assert joint.reached_genes == {"TGT", "TWO", "BIG", "ONE", "SUP_1", "INH", "MIX"}
    assert set(joint.genes) == joint.reached_genes
    assert joint.genes["TGT"].elements == ["E1"]  # the experimental twin was removed before gating
    for gene in ("TWO", "BIG", "ONE", "MIX", "TGT"):
        fr = joint.genes[gene].forms[sa.FORM]
        assert fr["ok"] and fr["harness_miss"] == [], gene
        assert fr["max_abs_run_minus_closed_form_log2"] <= sa.CLOSED_FORM_TOLERANCE
    assert joint.genes["BIG"].forms[sa.FORM]["summed_drop"] == pytest.approx(1.5)
    assert joint.genes["SUP_1"].forms[sa.FORM] == {"ok": False, "reason": "no_activating_rule"}
    assert joint.genes["INH"].forms[sa.FORM]["reason"] == "inhibitor_strength_above_1"
    assert joint.counts["runs"] == 5


def test_a_k562_hold_is_the_elements_own_observation_whatever_the_others_states(joint):
    two = joint.genes["TWO"]
    assert two.k562 == {"E5": "open", "E6": "closed"}
    change = two.forms[sa.FORM]["change"]
    assert set(change) == {"E5"}  # E6 is closed in K562: no hold run
    # under the mean and competition E6's hold moved E5's K562 change off its observation; here it does not
    assert change["E5"] == pytest.approx(0.4, abs=sa.CLOSED_FORM_TOLERANCE)
    mix = joint.genes["MIX"].forms[sa.FORM]["change"]
    assert mix["E10"] == pytest.approx(1.5, abs=sa.CLOSED_FORM_TOLERANCE)
    assert mix["E11"] == pytest.approx(-0.5, abs=sa.CLOSED_FORM_TOLERANCE)


# ---- the pairs through the multiplicative stages ---------------------------------------------------------


def test_a_pair_m2d_dropped_for_its_basal_rate_is_kept_here(traced):
    f, why = traced["big"]
    assert f["s4_bridge_kept"] == 0.0 and why["s4_bridge"] == "not_identifiable"
    assert f[f"{sa.S4}_kept"] == 1.0 and f[f"{sa.S4}_drop"] == pytest.approx(2.0)
    assert f[f"{sa.S5}_kept"] == 1.0 and f[f"{sa.S7}_kept"] == 1.0
    assert f[f"{sa.S7}_drop"] == pytest.approx(2.0, abs=sa.CLOSED_FORM_TOLERANCE)  # BY_CONSTRUCTION (3)
    f, _why = traced["two"]
    assert f[f"{sa.S4}_drop"] == pytest.approx(0.4) and f[f"{sa.S6}_kept"] == 1.0
    assert f[f"{sa.S7}_drop"] == pytest.approx(f[f"{sa.S5}_drop"], abs=sa.CLOSED_FORM_TOLERANCE)


def test_each_drop_names_its_reason(traced):
    f, why = traced["two6"]
    assert f[f"{sa.S4}_kept"] == 1.0 and f[f"{sa.S5}_kept"] == 0.0
    assert why[sa.S5] == "reader_closed" and why[sa.S7] == "held in K562"
    f, why = traced["sup"]
    assert f[f"{sa.S4}_kept"] == 0.0 and why[sa.S4] == "no_activating_rule"
    f, why = traced["inh"]
    assert f[f"{sa.S4}_kept"] == 0.0 and why[sa.S4] == "inhibitor_strength_above_1"
    f, why = traced["one"]
    assert f[f"{sa.S7}_drop"] == pytest.approx(1.0, abs=sa.CLOSED_FORM_TOLERANCE)
    assert why[sa.S6] == "no factor cluster names K562 at the element"
    f, why = traced["mix"]
    assert f[f"{sa.S4}_drop"] == pytest.approx(1.5)  # the inhibitor's 0 does not lower the maximum
    f, why = traced["none"]
    assert not any(f[f"{s}_kept"] for s in sa.MULT_STAGES) and sa.S4 not in why


def test_the_columns_never_read_the_label(look, joint):
    p, q = copy.deepcopy(PAIRS["big"]), copy.deepcopy(PAIRS["big"])
    q.regulated = not p.regulated
    assert sa.trace(p, look, joint) == sa.trace(q, look, joint)


# ---- the harness and the census --------------------------------------------------------------------------


def test_the_harness_passes_a_sound_trace_and_names_each_defect(look, joint, traced):
    rows = []
    for k, p in PAIRS.items():
        p.features.update(traced[k][0])
        rows.append(p)
    assert st.harness({"chrT": look}, rows) == []
    assert sa.harness({"chrT": look}, {"chrT": joint}, rows) == []
    bad = copy.deepcopy(joint)
    bad.genes["TWO"].forms[sa.FORM]["harness_miss"] = ["K562 hold E5: run 1, closed form 2"]
    bad.genes["TWO"].k562["E5"] = "closed"
    del bad.genes["ONE"]
    problems = sa.harness({"chrT": look}, {"chrT": bad}, [])
    assert any("not read from the gated program" in x for x in problems)
    assert any("K562 state closed against open" in x for x in problems)
    assert any("K562 hold E5" in x for x in problems)
    stray = pair(1, 2, "X")
    stray.features.update({f"{s}_{c}": 0.0 for s in sa.STAGES for c in ("kept", "drop")})
    stray.features[f"{sa.S7}_kept"] = 1.0
    assert sa.harness({}, {}, [stray]) == [f"{sa.S7} keeps 1 pairs that {sa.S5} dropped"]


def test_the_census_counts_the_identifiable_genes_and_the_concave_bound(joint):
    m2d = {
        "descriptives_not_verdicts": {"per_form": {"mean": {"runnable": 3}, "competition": {"runnable": 1}}}
    }
    c = sa.census({"chrT": joint}, m2d)
    assert c["reached_genes"] == 7 and c["identifiable"] == 5
    assert c["not_identifiable_by_reason"] == {"no_activating_rule": 1, "inhibitor_strength_above_1": 1}
    assert c["admitted_by_combine2_classify"] == 6  # SUP_1 is excluded by combine2.classify
    # TWO (0.37), ONE (0.5), MIX (0.65), TGT (0.29), INH (0.5) sum below 1; BIG (1.5) does not
    assert c["concave_bound_summed_drop_below_1"] == 5
    assert c["sentence"] == (
        "under the multiplicative composition 5 of the 7 genes the traced pairs reach are identifiable (2 are"
        " not: 1 inhibitor_strength_above_1, 1 no_activating_rule); M2d ran 3 of them under the runtime's"
        " mean and 1 under promoter competition, and a concave composition with a positive basal rate,"
        " whatever its parameters, could fit at most the 5 whose activating single-removal drops sum to"
        " less than 1"
    )


# ---- the readings ----------------------------------------------------------------------------------------


def _row(outcome: str) -> dict:
    return {"outcome": outcome, "reading": f"<{outcome}>"}


def test_the_chain_and_the_bridge_answer_follow_signaltraces_walk():
    stages = {
        "s2_target": _row("WEAKENED"),
        "s3_compiled_value": _row("WEAKENED"),
        sa.S4: _row("WEAKENED"),
        sa.S5: _row("UNDECIDED"),
        sa.S7: _row("UNDECIDED"),
    }
    kind, text, at = sa.chain_verdict(stages)
    assert (kind, at) == ("UNDECIDED_AT", sa.S5)
    assert text.endswith("<UNDECIDED>; it was first weakened at s2_target (target choice)")
    b, btext = sa.bridge_verdict(stages)
    assert b == "SURVIVES" and btext.startswith(
        "with the genes fitted jointly under the multiplicative composition (Zhou et al. 2024), the measured"
        " signal survives the bridge"
    )
    stages[sa.S4] = _row("NONE_KEPT")
    assert sa.bridge_verdict(stages)[0] == "LOST"
    kind, _text, at = sa.chain_verdict(stages)
    assert (kind, at) == ("LOST_AT", sa.S4)


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


def test_score_reads_the_stages_against_the_anchor_and_the_comparators():
    fit, pop = _population(seed=1), _population(seed=2)
    for p in fit + pop:
        f = p.features
        for s in sa.STAGES:
            f[f"{s}_kept"], f[f"{s}_drop"] = 0.0, 0.0
        for s in ("s1_anchor", "s2_target", "s3_compiled_value", sa.S4):
            f[f"{s}_kept"], f[f"{s}_drop"] = f["top_target"], f["deletion_drop"]
    anchor = st.anchor_reading(fit, pop)
    rows, scores, weights = sj.comparators(fit, pop, anchor)
    out = sa.score(fit, pop, anchor, rows, scores, weights)
    s = out["stages"]
    assert s[sa.S4]["gain_over_activity_distance"] == s["s2_target"]["gain_over_activity_distance"]
    assert s[sa.S4]["against_the_stage_it_follows"]["gain"] == 0.0
    assert s[sa.S5]["outcome"] == "NONE_KEPT" and out["chain"] == "LOST_AT"
    assert out["bridge"] == sa.BRIDGE_OF[s[sa.S4]["outcome"]]
    assert out["reading"].startswith("with the genes fitted jointly under the multiplicative composition")
    assert out["chain_reading"] in out["reading"]


# ---- the committed inputs this lane binds to ------------------------------------------------------------


def test_m2ds_committed_census_is_what_the_census_quotes():
    import tracked_paths as tp

    path = Path("data/results/s4joint_result.json")
    tp.must_be_present(path)
    per = json.loads(path.read_text())["descriptives_not_verdicts"]["per_form"]
    assert {f: per[f]["reached_genes"] for f in ("mean", "competition")} == {"mean": 210, "competition": 210}
    assert {f: per[f]["runnable"] for f in ("mean", "competition")} == {"mean": 10, "competition": 2}


def test_the_committed_registration_is_the_module_it_binds_and_reuses_signaltraces_population():
    import tracked_paths as tp

    reg_path = Path("data/results/saturate_registration.json")
    st_path = Path("data/results/signaltrace_registration.json")
    tp.must_be_present(reg_path, st_path)
    reg = json.loads(reg_path.read_text())
    st_reg = json.loads(st_path.read_text())
    assert {f: reg[f] for f in sa.BOUND} == sa.registered_fields()
    assert reg["population_record"]["sha256"] == st_reg["population_sha256"]
    assert reg["population_record"]["pairs"] == len(st_reg["population_ids"]) == 1918
    assert reg["comparator_committed_verbatim"]["s4_bridge"].startswith("none kept")
    assert {f: v["runnable"] for f, v in reg["m2d_committed_verbatim"]["census"].items()} == {
        "mean": 10,
        "competition": 2,
    }
    assert reg["result_manifest"]["complete"] is True
    assert reg["result_manifest"]["code_cleanliness"]["own_code_is_committed"] is True
    assert reg["alphagenome_requests"] == 0
