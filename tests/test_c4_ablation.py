# SPDX-License-Identifier: AGPL-3.0-or-later
"""The AlphaGenome ablation (item 13 C4): the classification rule on synthetic links, programs and
results (no data/knowledge needed)."""

from __future__ import annotations

import json

from genomeos.attribution import ablation as ab
from genomeos.attribution import measured as ms


def _pair(start, gene, outcome, cell="K562", split=ms.TRAINING, end=None):
    significant = outcome in (ms.DECREASE, ms.INCREASE)
    effect = -0.3 if outcome == ms.DECREASE else 0.3 if outcome == ms.INCREASE else -0.01
    power = 0.9 if outcome == ms.NULL_INFORMATIVE else 0.1
    return ms.CrispriPair(
        chrom="chr1",
        start=start,
        end=end or start + 300,
        gene=gene,
        cell=cell,
        dataset="S",
        reference="",
        regulated=outcome == ms.DECREASE,
        significant=significant,
        effect_size=effect,
        p_adjusted=0.01,
        split=split,
        power_at_effect_size_20=power,
    )


def test_a_link_survives_only_on_the_same_gene_in_the_stated_direction():
    assert ab.link_class([_pair(0, "G", ms.DECREASE)], "G", "activates")[0] == ab.SURVIVES
    assert ab.link_class([_pair(0, "G", ms.INCREASE)], "G", "inhibits")[0] == ab.SURVIVES
    # the opposite direction and a well-powered null contradict
    assert ab.link_class([_pair(0, "G", ms.INCREASE)], "G", "activates")[1] == ab.CONTRADICTED
    assert ab.link_class([_pair(0, "G", ms.NULL_INFORMATIVE)], "G", "activates")[1] == ab.CONTRADICTED
    # an underpowered null decides nothing
    assert ab.link_class([_pair(0, "G", ms.NULL_INCONCLUSIVE)], "G", "activates")[1] == ab.INCONCLUSIVE
    # another gene regulated is not this link measured
    cls, reason, detail = ab.link_class([_pair(0, "OTHER", ms.DECREASE)], "G", "activates")
    assert (cls, reason) == (ab.DOES_NOT_SURVIVE, ab.NEVER_MEASURED) and "not tested" in detail
    assert ab.link_class([], "G", "activates")[2] == "element never screened"


def test_a_rule_survives_only_in_its_own_cell():
    pairs = [_pair(0, "G", ms.DECREASE, cell="K562")]
    assert ab.link_class(pairs, "G", "activates", "K562")[0] == ab.SURVIVES
    assert ab.link_class(pairs, "G", "activates", "HepG2")[1] == ab.NEVER_MEASURED


def test_the_overlap_rule_is_the_measured_layers():
    idx = ab.CrispriIndex([_pair(1000, "G", ms.DECREASE, end=1100)])
    assert idx.of("chr1", 1000, 1100) and not idx.of("chr1", 900, 1400)  # 100 of 500 bp: below 0.5


def test_axis_values_sort_into_the_three_classes():
    link = (ab.SURVIVES, None)
    assert ab.axis_class("origin", "unique", False, link) == (ab.NEVER_MODEL_DEPENDENT, None)
    assert ab.axis_class("molecular_role", "enhancer_like", False, link)[0] == ab.NEVER_MODEL_DEPENDENT
    assert ab.axis_class("molecular_role", ab.REPRESSION_GROUP, False, link) == link
    assert ab.axis_class("activity", "activates_target", False, link) == link
    about = (ab.DOES_NOT_SURVIVE, ab.ABOUT_THE_MODEL)
    assert ab.axis_class("evidence_status", "predicted_model", False, link) == about
    assert ab.axis_class("evidence_status", "conflicting", True, link) == about
    # a measured twin's activity comes from the assays
    assert ab.axis_class("activity", "activates_target", True, link)[0] == ab.NEVER_MODEL_DEPENDENT


PROGRAM = """element E1 {
  class: enhancer
  locus: chr1:1000-1300
  targets: G1
  origin: unique
  molecular_role: enhancer_like
  activity: activates_target
  target_relation: predicted_deletion_target
  evidence_status: registry_biochemical, predicted_model, selection_not_measured
}
rule E1 activates G1 { strength: 0.4; when: cell_type = HepG2; evidence: predicted "x" effect -0.4 }
element E2 {
  class: enhancer
  locus: chr1:9000-9300
  targets: G2
  origin: unique
  molecular_role: enhancer_like, silencer|insulator_like|competing_promoter|unknown
  activity: represses_target
  target_relation: predicted_deletion_target
  evidence_status: predicted_model
}
rule E2 inhibits G2 { strength: 0.2; when: cell_type = K562; evidence: predicted "x" effect 0.2 }
element E1_measured {
  class: enhancer
  locus: chr1:1000-1300
  activity: activates_target
  target_relation: measured_perturbation_target
  evidence_status: measured
}
rule E1_measured activates G1 { strength: 0.3; when: cell_type = K562; evidence: experimental "y" }
"""
PROGRAM += (
    'rule E1_measured activates G7 { strength: 0.3; when: cell_type = K562; evidence: experimental "y, '
    + ms.HELDOUT_MARK
    + '" }\n'
)


def test_the_program_census_classifies_elements_rules_and_axes(tmp_path):
    p = tmp_path / "noncoding_chr1.bio"
    p.write_text(PROGRAM)
    idx = ab.CrispriIndex([_pair(1000, "G1", ms.DECREASE), _pair(9000, "G2", ms.NULL_INFORMATIVE)])
    c = ab.program_census([p], idx)
    assert c["elements"][(ab.SURVIVES, None)] == 1
    assert c["elements"][(ab.DOES_NOT_SURVIVE, ab.CONTRADICTED)] == 1
    # E1's rule is gated on HepG2 but was measured in K562: the element's link survives, its rule does not
    assert c["rules"][(ab.DOES_NOT_SURVIVE, ab.NEVER_MEASURED)] == 1
    assert c["rules"][(ab.DOES_NOT_SURVIVE, ab.CONTRADICTED)] == 1
    assert c["rules"][(ab.NEVER_MODEL_DEPENDENT, None)] == 2
    assert c["experimental_rules"]["heldout_marked"] == 1
    pred = c["axes"]["predicted"]
    assert pred["evidence_status"][("predicted_model", ab.DOES_NOT_SURVIVE, ab.ABOUT_THE_MODEL)] == 2
    assert pred["molecular_role"][(ab.REPRESSION_GROUP, ab.DOES_NOT_SURVIVE, ab.CONTRADICTED)] == 1
    assert pred["origin"][("unique", ab.NEVER_MODEL_DEPENDENT, None)] == 2
    assert c["axes"]["measured"]["activity"][("activates_target", ab.NEVER_MODEL_DEPENDENT, None)] == 1


def test_placement_counts_measured_links_the_model_compiled_no_element_for():
    pairs = [
        _pair(1000, "G1", ms.DECREASE),
        _pair(50_000, "G2", ms.DECREASE),
        _pair(80_000, "G3", ms.DECREASE),
    ]
    out = ab.placement(pairs, {"chr1": [(1000, 1300)]}, {"chr1": [(50_000, 50_300)]})
    assert out["measured_links"] == 3
    assert out["placed_training"] == 1
    assert out["not_placed_registry_would_training"] == 1
    assert out["not_placed_no_registry_element_training"] == 1


def _write(d, name, payload):
    (d / f"{name}.json").write_text(json.dumps(payload))


def test_the_headline_census_applies_the_registered_readings(tmp_path):
    mdl = {"result_manifest": {"inputs": [{"path": "data/knowledge/alphagenome/all_elements/chr1.json"}]}}
    uniform = {"excess_points": 5.0, "bootstrap_chromosomes": {"ci95": [1.0, 9.0]}}
    _write(tmp_path, "manifest_headlines", {"rebuilt": [{"result": n} for n in ab.HEADLINES]})
    _write(
        tmp_path,
        "node_containment_audit",
        {**mdl, "modelled": {"controls": {"uniform": {"excess_points": 2.9}}}},
    )
    _write(tmp_path, "node_containment_measured", {**mdl, "measured": {"controls": {"uniform": uniform}}})
    _write(tmp_path, "constrained_unknown_targets", mdl)
    _write(tmp_path, "clause2_measured_arm", {"primary_reading": {"outcome": "cannot_decide"}})
    _write(
        tmp_path, "therapeutic_benchmark", {"result_manifest": {"inputs": [{"path": "data/knowledge/x.tsv"}]}}
    )
    targets = {"result_manifest": {"inputs": [{"path": "data/results/enhancer_targets_all_chr1.json"}]}}
    _write(tmp_path, "unknown_coverage", targets)
    _write(tmp_path, "crispri_published", mdl)
    rows = {r["result"]: r for r in ab.headline_census(tmp_path)}
    assert rows["node_containment_audit"]["cls"] == ab.SURVIVES
    assert rows["node_containment_audit"]["arm"]["value"] == 5.0
    assert rows["node_containment_measured"]["cls"] == ab.NEVER_MODEL_DEPENDENT
    cut = rows["constrained_unknown_targets"]
    assert (cut["cls"], cut["reason"]) == (ab.DOES_NOT_SURVIVE, ab.INCONCLUSIVE)
    assert rows["therapeutic_benchmark"]["cls"] == ab.NEVER_MODEL_DEPENDENT
    assert rows["unknown_coverage"]["cls"] == ab.NEVER_MODEL_DEPENDENT
    assert rows["crispri_published"]["reason"] == ab.ABOUT_THE_MODEL


def test_an_arm_whose_interval_crosses_zero_does_not_carry_the_conclusion(tmp_path):
    mdl = {"result_manifest": {"inputs": [{"path": "alphagenome/x"}]}}
    uniform = {"excess_points": 1.0, "bootstrap_chromosomes": {"ci95": [-1.0, 3.0]}}
    _write(tmp_path, "manifest_headlines", {"rebuilt": [{"result": "node_containment_audit"}]})
    _write(
        tmp_path,
        "node_containment_audit",
        {**mdl, "modelled": {"controls": {"uniform": {"excess_points": 2.9}}}},
    )
    _write(tmp_path, "node_containment_measured", {"measured": {"controls": {"uniform": uniform}}})
    (row,) = ab.headline_census(tmp_path)
    assert (row["cls"], row["reason"]) == (ab.DOES_NOT_SURVIVE, ab.INCONCLUSIVE)


def test_a_model_free_registration_contradicted_by_the_manifest_is_not_classified(tmp_path):
    bad = {"result_manifest": {"inputs": [{"path": "data/knowledge/alphagenome/e.json"}]}}
    _write(tmp_path, "manifest_headlines", {"rebuilt": [{"result": "therapeutic_benchmark"}]})
    _write(tmp_path, "therapeutic_benchmark", bad)
    (row,) = ab.headline_census(tmp_path)
    assert row["cls"] == "registration_contradicted"
