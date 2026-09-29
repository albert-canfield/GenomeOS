# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review R3 (2026-09-28): annotation to dynamics, registered in 06b7557 before the build.

Negatives first: what the runtime must say when a model is not resolved, what it must not say when a
caller deliberately reads zero, and that another citation of one mechanism is not another parameter.
Then the end-to-end fixture: one compiled deletion observation reproduced by the simulation it maps to.
"""

from __future__ import annotations

import json
import math
import warnings

import pytest

from genomeos.attribution import bridge
from genomeos.attribution.compile import compile_chromosome
from genomeos.attribution.measured import Layer
from genomeos.lang.parser import parse
from genomeos.runtime.grn import (
    EXPLICIT_ZERO_SUFFIX,
    NetworkRuntime,
    UnresolvedModel,
    UnresolvedModelWarning,
)

K562 = {"cell_type": "K562"}


def _steady(module, context, clamp, hours=30.0):
    rt = NetworkRuntime(module, context=context, strict=True)
    return rt.run(hours=hours, dt=0.01, clamp=clamp, record_every=100).final()


# ---- negatives ------------------------------------------------------------------------------------


def test_missing_regulator_state_is_an_explicit_unresolved_model():
    m = parse(
        "element X { class: enhancer }\ngene G { basal: 1; max: 10 }\nrule X activates G { strength: 0.5 }\n",
        "t",
    )
    with pytest.raises(UnresolvedModel) as err:
        NetworkRuntime(m, strict=True).run(hours=1.0)
    assert [(i.reason, i.subject) for i in err.value.issues] == [("regulator_state_missing", "X")]
    with pytest.warns(UnresolvedModelWarning):
        traj = NetworkRuntime(m).run(hours=1.0)
    assert [i.reason for i in traj.unresolved] == ["regulator_state_missing"]
    # the value integrated is unchanged: the diagnostic reports, it does not repair
    assert traj.final()["G.mRNA"] == pytest.approx(1 - math.exp(-1.0), rel=1e-3)


def test_a_clamped_regulator_is_resolved():
    m = parse(
        "element X { class: enhancer }\ngene G { basal: 1; max: 10 }\nrule X activates G { strength: 0.5 }\n",
        "t",
    )
    traj = NetworkRuntime(m, strict=True).run(hours=1.0, clamp={"X": 1.0})
    assert traj.unresolved == []


def test_a_declared_zero_is_not_unresolved_and_reads_zero():
    """network_experiment's edge knockout renames the source '<id>@zero': deliberate, so not unresolved."""
    m = parse("gene X { basal: 1 }\ngene G { basal: 1; max: 10 }\nrule X activates G { strength: 1 }\n", "t")
    m.rules[0].source = f"X{EXPLICIT_ZERO_SUFFIX}"
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        traj = NetworkRuntime(m, strict=True).run(hours=20.0)
    assert traj.unresolved == []
    assert traj.final()["G.mRNA"] == pytest.approx(1.0, rel=1e-3)  # basal only: the source read 0


def test_a_regulated_gene_without_max_rate_is_reported():
    m = parse(
        "element X { class: enhancer }\ngene G { basal: 1 }\nrule X activates G { strength: 0.5 }\n", "t"
    )
    with pytest.raises(UnresolvedModel) as err:
        NetworkRuntime(m, strict=True).run(hours=1.0, clamp={"X": 1.0})
    assert [(i.reason, i.subject) for i in err.value.issues] == [("gene_parameter_missing", "G")]


def _compiled(tmp_path, lfc=-1.0, extra=""):
    block = {
        "start": 0,
        "end": 1000,
        "length": 1000,
        "class": "regulatory",
        "phylop": None,
        "elements": None,
        "guess": {"tier": "regulatory", "label": "regulatory elements", "confidence": 0.5},
    }
    (tmp_path / "budget_chrT.json").write_text(
        json.dumps({"chrom": "chrT", "unknown_bp": 1000, "constrained_fraction": 0.0, "blocks": [block]})
    )
    el = {
        "id": "E1",
        "start": 100,
        "end": 200,
        "predicted_coding": {
            "gene": "GENE1",
            "action": "activates" if lfc < 0 else "represses",
            "log2_fold_change": lfc,
            "tissue": "K562",
            "strength": "strong",
            "confidence": abs(lfc),
        },
    }
    (tmp_path / "enhancer_targets_chrT.json").write_text(json.dumps({"elements": [el]}))
    text = compile_chromosome("chrT", tmp_path, layer=Layer(chrom="chrT")) + extra
    return parse(text, "chrT")


def test_a_compiled_program_as_it_stands_is_annotation_and_says_so(tmp_path):
    m = _compiled(tmp_path)
    issues = NetworkRuntime(m, context=K562).diagnose()
    assert {(i.reason, i.subject) for i in issues} == {
        ("regulator_state_missing", "E1"),
        ("gene_parameter_missing", "GENE1"),
    }


def test_a_censored_observation_is_unresolved():
    m = parse(
        "element E2 { class: enhancer }\ngene G { basal: 1; max: 10 }\n"
        'rule E2 activates G { strength: 1.0; when: cell_type = K562; evidence: predicted "AlphaGenome'
        ' deletion, K562" }\n',
        "t",
    )
    b = bridge.parameterize(m, K562)
    assert [(i.reason, i.subject) for i in b.unresolved] == [("observation_censored", "E2")]
    assert b.strengths == {}


def test_a_gene_without_declared_parameters_is_unresolved(tmp_path):
    b = bridge.parameterize(_compiled(tmp_path), K562)  # the compiled stub declares neither
    assert [(i.reason, i.subject) for i in b.unresolved] == [("gene_parameter_missing", "GENE1")]


def test_two_mechanisms_on_one_gene_are_not_identifiable(tmp_path):
    extra = (
        "element E9 { class: enhancer }\n"
        'rule E9 activates GENE1 { strength: 0.3; when: cell_type = K562; evidence: predicted "AlphaGenome'
        ' deletion, K562" effect -0.3 log2 fold change on deletion, probability unavailable }\n'
    )
    b = bridge.parameterize(_compiled(tmp_path, extra=extra), K562, genes=GENES)
    assert [i.reason for i in b.unresolved] == ["not_identifiable"] and b.strengths == {}


def test_a_response_the_declared_parameters_cannot_reach_is_reported_not_clipped(tmp_path):
    b = bridge.parameterize(_compiled(tmp_path, lfc=-6.0), K562, genes=GENES)  # 64-fold loss, V = 10 b
    assert [i.reason for i in b.unresolved] == ["response_out_of_range"] and b.strengths == {}


def test_under_the_mean_rule_adding_an_activator_can_lower_expression():
    """A declared model assumption (lane-sign, b1f3405): not a bridge decision, pinned so it stays visible."""
    head = "element A { class: enhancer }\nelement B { class: enhancer }\ngene G { basal: 0; max: 10 }\n"
    one = parse(head + "rule A activates G { strength: 1 }\n", "t")
    two = parse(head + "rule A activates G { strength: 1 }\nrule B activates G { strength: 1 }\n", "t")
    g1 = _steady(one, {}, {"A": 5.0})["G.mRNA"]
    g2 = _steady(two, {}, {"A": 5.0, "B": 0.0})["G.mRNA"]
    assert g2 < 0.6 * g1 and bridge.ACTIVATOR_COMBINATION == "mean"


# ---- citations: one mechanism, one parameter ------------------------------------------------------

MEASURED = (
    "element E1_measured {{ class: enhancer; evidence: experimental"
    ' "CRISPRi enhancer-gene screens" links GENE1 in K562 from {n} pairs }}\n'
    "rule E1_measured activates GENE1 {{ strength: {f}; when: cell_type = K562; evidence: experimental"
    ' "CRISPRi enhancer-gene screens, silenced in K562"; confidence: 0.9 }}\n'
)
DUPLICATE = (
    'rule E1 activates GENE1 { strength: 1.0; when: cell_type = K562; evidence: predicted "another'
    ' deletion study, K562" effect -1 log2 fold change on deletion, probability unavailable }\n'
)
GENES = {"GENE1": {"basal_rate": 1.0, "max_rate": 10.0}}


def test_another_citation_of_the_same_mechanism_does_not_change_its_strength(tmp_path):
    alone = bridge.parameterize(_compiled(tmp_path), K562, genes=GENES)
    ((key, s),) = alone.strengths.items()
    dup = bridge.parameterize(_compiled(tmp_path, extra=DUPLICATE), K562, genes=GENES)
    assert dup.strengths == {key: s} and not dup.unresolved
    only = parse(MEASURED.format(f=0.3, n=1) + "gene GENE1 { symbol: GENE1 }\n", "m")
    measured = bridge.parameterize(only, K562, genes=GENES)
    both = bridge.parameterize(_compiled(tmp_path, extra=MEASURED.format(f=0.3, n=1)), K562, genes=GENES)
    assert both.strengths == measured.strengths  # the predicted rule beside it adds nothing
    assert [r.source for r in both.superseded] == ["E1"]
    assert len(both.module.rules) == 1  # one rule per mechanism enters the simulation


def test_conflicting_citations_of_one_kind_are_reported_not_averaged(tmp_path):
    other = DUPLICATE.replace("effect -1 log2", "effect -2 log2")
    b = bridge.parameterize(_compiled(tmp_path, extra=other), K562, genes=GENES)
    assert [i.reason for i in b.unresolved] == ["conflicting_observations"] and b.strengths == {}


# ---- the end-to-end fixture -----------------------------------------------------------------------


@pytest.mark.parametrize("lfc", [-1.0, -0.4, 0.7])
def test_a_compiled_deletion_is_reproduced_by_the_simulation_it_maps_to(tmp_path, lfc):
    """Acceptance (1): compile, parse, parameterise, run intact and removed; the ratio is the observation."""
    b = bridge.parameterize(_compiled(tmp_path, lfc=lfc), K562, genes=GENES)
    assert not b.unresolved and b.required_state == {"E1": bridge.INTACT}
    intact = _steady(b.module, K562, {"E1": bridge.INTACT})["GENE1.mRNA"]
    removed = _steady(b.module, K562, {"E1": bridge.REMOVED})["GENE1.mRNA"]
    assert removed / intact == pytest.approx(2.0**lfc, rel=0.01)
    # and the strength the compiler wrote, |lfc|, is not the parameter that reproduces it
    assert not math.isclose(next(iter(b.strengths.values())), abs(lfc), rel_tol=0.01)


def test_a_measured_crispri_effect_is_reproduced_as_one_plus_f():
    m = parse(MEASURED.format(f=0.3, n=1) + "gene GENE1 { symbol: GENE1 }\n", "m")
    b = bridge.parameterize(m, K562, genes=GENES)
    intact = _steady(b.module, K562, {"E1": 1.0})["GENE1.mRNA"]
    removed = _steady(b.module, K562, {"E1": 0.0})["GENE1.mRNA"]
    assert removed / intact == pytest.approx(0.7, rel=0.01)


def test_a_crispri_strength_taken_as_the_largest_of_several_pairs_is_not_fitted():
    """The amendment of b986238 (census A8): a maximum over pairs is reported, never fitted to."""
    m = parse(MEASURED.format(f=0.3, n=3) + "gene GENE1 { symbol: GENE1 }\n", "m")
    b = bridge.parameterize(m, K562, genes=GENES)
    assert [(i.reason, i.subject) for i in b.unresolved] == [("observation_summarised", "E1_measured")]
    assert "largest of 3 pairs" in b.unresolved[0].detail and b.strengths == {}
