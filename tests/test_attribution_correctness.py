# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item 12 S4: what a correct attribution means. The assay-by-axis table, the scorer that keeps target
accuracy, role accuracy and coverage apart, the observation model's inadequacy as an outcome, and
holdout's split discipline, on synthetic units only (no data/knowledge needed)."""

from __future__ import annotations

import json
from collections import Counter, defaultdict

import pytest

from genomeos.attribution import correctness as co
from genomeos.attribution import holdout as ho
from genomeos.attribution import measured as ms

E = ("chr1", 10_000, 10_400)  # the element every synthetic claim is about


def claim(axis: str, value: str, gene: str = "", cell: str = "") -> co.Claim:
    if axis == co.TARGET:
        gene = gene or value
    return co.Claim("EH1", *E, axis, value, gene, cell)


def obs(kind: str, gene: str = "G1", cell: str = "K562", **kw) -> co.Observation:
    return co.Observation(
        kind, kw.pop("source", f"src:{kind}"), gene if kind in co.PAIR_KINDS else "", cell, **kw
    )


def unit(source: str, outcome: str, gene: str = "G1", cell: str = "K562", **kw) -> ho.Unit:
    start, end = kw.pop("start", E[1]), kw.pop("end", E[2])
    return ho.Unit(
        source=source, chrom=E[0], start=start, end=end, outcome=outcome, gene=gene, cell=cell, **kw
    )


def labels(reads=frozenset({ho.MODEL}), built_without=None) -> ho.Labels:
    return ho.Labels("synthetic", lambda u, e: None, frozenset(reads), built_without=built_without)


# --- the table ------------------------------------------------------------------------------------
def test_every_observation_kind_has_one_cell_per_axis():
    assert set(co.TABLE) == set(co.OBSERVATIONS)
    for kind in co.OBSERVATIONS:
        assert set(co.TABLE[kind]) == set(co.AXES), kind
        assert co.OBSERVATION_TEXT[kind]
        assert (kind in co.LOADED) != (kind in co.NOT_LOADED), kind
    assert set(co.CRISPRI_KIND) == set(ms.OUTCOMES)  # R2's five outcomes, each a row


def test_the_reviews_cannot_cells():
    # a reporter tile can establish activity, never a target
    assert co.TABLE["lentimpra_active"][co.TARGET].kind == co.CANNOT
    assert co.reading("lentimpra_active", co.ACTIVITY, "active_in_reporter") == co.ESTABLISHES
    # reporter activity never substitutes for regulation in place
    for v in co.DIRECTION:
        assert co.reading("lentimpra_active", co.ACTIVITY, v) == co.CANNOT
        assert co.reading("vista_positive", co.ACTIVITY, v) == co.CANNOT
    # conservation can suggest selection, never a role (nor anything else)
    for a in co.AXES:
        assert co.TABLE["conservation_constrained"][a].kind == co.CANNOT
        assert co.TABLE["alphagenome_deletion_prediction"][a].kind == co.CANNOT
        assert co.TABLE["crispri_null_underpowered"][a].kind == co.CANNOT
    # not naming a gene is not absence of function: a null says nothing about role or origin
    assert co.TABLE["crispri_null_well_powered"][co.ROLE].kind == co.CANNOT
    assert co.TABLE["crispri_null_well_powered"][co.ORIGIN].kind == co.CANNOT
    # an increase is not a silencer
    assert co.reading("crispri_increase", co.ROLE, "silencer") == co.SUGGESTS
    # contact and association never establish a target
    for k in ("chromatin_contact_present", "gtex_associated", "allelic_imbalance"):
        assert co.TABLE[k][co.TARGET].kind == co.SUGGESTS
    # nothing loaded judges origin or molecular role
    for k in co.LOADED:
        assert co.TABLE[k][co.ORIGIN].kind != co.JUDGES
        assert co.TABLE[k][co.ROLE].kind != co.JUDGES


def _values(axis: str) -> list[str]:
    if axis == co.ORIGIN:
        return list(co.ORIGIN_VALUES)
    if axis == co.ROLE:
        return [v for v in co.R7_AXES[co.ROLE] if v != "unknown"]
    if axis == co.ACTIVITY:
        return [v for v in co.R7_AXES[co.ACTIVITY] if v != "unknown"]
    return ["G1"] if axis == co.TARGET else ["K562"]


@pytest.mark.parametrize("kind", co.OBSERVATIONS)
def test_an_assay_never_judges_an_axis_the_table_says_it_cannot(kind):
    """Observations of one kind, in the claim's scope (same gene, same cell, same value), several of
    them: wherever the table says cannot or suggests, the claim stays not judged."""
    for axis in co.AXES:
        for value in _values(axis):
            r = co.reading(kind, axis, value)
            if r in (co.ESTABLISHES, co.REFUTES, co.READ_INADEQUATE):
                continue
            gene = "G1" if axis in (co.TARGET, co.CONTEXT, co.ACTIVITY) else ""
            c = claim(axis, value, gene=gene, cell="K562")
            seen = [obs(kind, value=value), obs(kind, value=value, source="another"), obs(kind, cell="HepG2")]
            v = co.verdict_of(c, seen)
            assert v.verdict == co.NOT_JUDGED, (kind, axis, value, v)
            assert v.reason in (co.SUGGESTS_ONLY, co.CANNOT_ONLY, co.OUTSIDE_SCOPE, co.CONDITION_UNMET)
            if axis == co.TARGET and co.TABLE[kind][axis].kind == co.CANNOT and kind not in co.SCREEN_KINDS:
                assert v.reason == co.CANNOT_ONLY


def test_only_cells_that_judge_decide_among_mixed_observations():
    seen = [obs(k) for k in co.OBSERVATIONS]
    for axis in co.AXES:
        for value in _values(axis):
            v = co.verdict_of(claim(axis, value, gene="G1", cell="K562"), seen)
            for o in v.deciding:
                assert co.reading(o.kind, axis, value) in (co.ESTABLISHES, co.REFUTES, co.READ_INADEQUATE)


# --- the axis rules -----------------------------------------------------------------------------
def test_response_elsewhere_and_null_in_the_stated_cell_is_one_context_error_not_two():
    seen = [obs("crispri_decrease", cell="HepG2"), obs("crispri_null_well_powered", cell="K562")]
    t = co.verdict_of(claim(co.TARGET, "G1", cell="K562"), seen)
    c = co.verdict_of(claim(co.CONTEXT, "K562", gene="G1", cell="K562"), seen)
    # S4's direction rule as registered, v1 (versioned since 2026-09-29; v2 has its own tests)
    a = co.verdict_of(claim(co.ACTIVITY, "activates_target", gene="G1", cell="K562"), seen, rule=co.RULE_V1)
    assert (t.verdict, c.verdict, a.verdict) == (co.CORRECT, co.INCORRECT, co.CORRECT)
    assert t.refutable and c.refutable


def test_a_null_in_the_stated_cell_refutes_the_target_and_leaves_role_and_context_unjudged():
    seen = [obs("crispri_null_well_powered", cell="K562")]
    t = co.verdict_of(claim(co.TARGET, "G1", cell="K562"), seen)
    assert t.verdict == co.INCORRECT
    assert co.OBSERVATION_MODEL_INADEQUATE in t.explanations
    assert any("redundant" in x for x in t.explanations)
    for c in (
        claim(co.CONTEXT, "K562", gene="G1", cell="K562"),
        claim(co.ACTIVITY, "activates_target", gene="G1", cell="K562"),
    ):
        v = co.verdict_of(c, seen)
        assert (v.verdict, v.reason) == (co.NOT_JUDGED, co.CONDITION_UNMET)
    for c in (claim(co.ROLE, "enhancer_like"), claim(co.ORIGIN, "unique")):
        assert co.verdict_of(c, seen).verdict == co.NOT_JUDGED


def test_a_null_only_in_another_cell_is_not_a_refutation():
    v = co.verdict_of(
        claim(co.TARGET, "G1", cell="placenta"), [obs("crispri_null_well_powered", cell="K562")]
    )
    assert (v.verdict, v.reason, v.detail) == (co.NOT_JUDGED, co.OUTSIDE_SCOPE, co.NULL_ELSEWHERE)
    assert not v.refutable


def test_a_claim_with_no_stated_cell_can_be_established_and_never_refuted():
    null = [obs("crispri_null_well_powered", cell="K562")]
    assert co.verdict_of(claim(co.TARGET, "G1", cell="unknown"), null).verdict == co.NOT_JUDGED
    resp = [obs("crispri_decrease", cell="K562")]
    assert co.verdict_of(claim(co.TARGET, "G1", cell="unknown"), resp).verdict == co.CORRECT


def test_direction_and_reasons():
    v1 = co.RULE_V1  # S4's direction rule as registered (versioned since 2026-09-29; v2 has its own tests)
    opposite = [obs("crispri_increase")]
    assert co.verdict_of(claim(co.TARGET, "G1", cell="K562"), opposite).verdict == co.CORRECT
    a = co.verdict_of(claim(co.ACTIVITY, "activates_target", gene="G1", cell="K562"), opposite, rule=v1)
    assert a.verdict == co.INCORRECT and co.OBSERVATION_MODEL_INADEQUATE in a.explanations
    split = [obs("crispri_decrease", cell="K562"), obs("crispri_increase", cell="HepG2")]
    v = co.verdict_of(claim(co.ACTIVITY, "activates_target", gene="G1", cell="WTC11"), split, rule=v1)
    assert v.verdict == co.UNRESOLVED and co.OBSERVATION_MODEL_INADEQUATE in v.explanations
    # in the stated cell, the stated cell decides
    assert (
        co.verdict_of(claim(co.ACTIVITY, "activates_target", gene="G1", cell="K562"), split, rule=v1).verdict
        == co.CORRECT
    )
    other_gene = [obs("crispri_decrease", gene="G2")]
    v = co.verdict_of(claim(co.TARGET, "G1", cell="K562"), other_gene)
    assert (v.reason, v.detail) == (co.OUTSIDE_SCOPE, co.GENE_NOT_TESTED)
    weak = [obs("crispri_null_underpowered")]
    assert co.verdict_of(claim(co.TARGET, "G1", cell="K562"), weak).reason == co.CANNOT_ONLY
    assert (
        co.verdict_of(claim(co.TARGET, "G1", cell="K562"), [obs("gtex_associated")]).reason
        == co.SUGGESTS_ONLY
    )
    tile = [obs("lentimpra_active")]
    assert co.verdict_of(claim(co.TARGET, "G1", cell="K562"), tile).reason == co.CANNOT_ONLY
    assert co.verdict_of(claim(co.ROLE, "enhancer_like"), tile).reason == co.SUGGESTS_ONLY
    assert co.verdict_of(claim(co.TARGET, "G1", cell="K562"), []).reason == co.NO_OBSERVATION


def test_unknown_and_unchosen_alternatives_are_not_claims():
    for value in ("unknown", "unassigned", "silencer|insulator_like|competing_promoter|unknown"):
        with pytest.raises(ValueError):
            co.verdict_of(claim(co.ROLE, value), [])


def test_value_kinds_judge_origin_and_role_only_when_they_state_a_value():
    rm = [obs("sequence_annotation", value="repeat_derived/LINE")]
    assert co.verdict_of(claim(co.ORIGIN, "repeat_derived/LINE"), rm).verdict == co.CORRECT
    assert co.verdict_of(claim(co.ORIGIN, "unique"), rm).verdict == co.INCORRECT
    assert co.verdict_of(claim(co.ORIGIN, "repeat_derived"), rm).verdict == co.CORRECT
    assert co.verdict_of(claim(co.ORIGIN, "repeat_derived/SINE"), rm).verdict == co.INCORRECT
    reg = [obs("registry_biochemical", value="enhancer_like,insulator_like")]
    assert co.verdict_of(claim(co.ROLE, "insulator_like"), reg).verdict == co.CORRECT
    assert co.verdict_of(claim(co.ROLE, "promoter_like"), reg).verdict == co.INCORRECT
    # nothing establishes a silencer
    assert co.verdict_of(claim(co.ROLE, "silencer"), reg).verdict == co.NOT_JUDGED


# --- the observation model is inadequate -----------------------------------------------------------
def test_the_observation_model_is_inadequate_is_a_verdict_and_survives_serialisation():
    tiles = [obs("lentimpra_tiles_conflict", cell="K562")]
    v1 = co.verdict_of(claim(co.ACTIVITY, "active_in_reporter", cell="K562"), tiles)
    both = [obs("crispri_decrease", cell="K562"), obs("crispri_increase", cell="K562", source="other")]
    v2 = co.verdict_of(claim(co.ACTIVITY, "represses_target", gene="G1", cell="K562"), both)
    for v in (v1, v2):
        assert v.verdict == co.MODEL_INADEQUATE
        assert v.explanations == (co.OBSERVATION_MODEL_INADEQUATE,)
        back = co.Verdict.from_dict(json.loads(json.dumps(v.to_dict())))
        assert back == v
    assert co.MODEL_INADEQUATE in co.VERDICTS and co.MODEL_INADEQUATE in co.JUDGED
    for verdict in (co.INCORRECT, co.UNRESOLVED, co.CORRECT):
        assert co.OBSERVATION_MODEL_INADEQUATE in co.EXPLANATIONS[verdict]


# --- the three quantities ------------------------------------------------------------------------
def _synthetic_units() -> dict[str, tuple[ho.Unit, ...]]:
    far = {"start": 900_000, "end": 900_400}
    return {
        "crispri:A": (
            unit("crispri:A", ms.DECREASE, split=ms.TRAINING, study="A"),
            unit("crispri:A", ms.NULL_INFORMATIVE, gene="G2", split=ms.TRAINING, study="A"),
            unit("crispri:A", ms.INCREASE, gene="G3", split=ms.HELDOUT, study="A", **far),
        ),
        "lentimpra:K562": (
            unit("lentimpra:K562", "active", gene="", value=1.5),
            unit("lentimpra:K562", "inactive", gene="", value=0.1, **far),
        ),
        "vista": (unit("vista", "positive", gene="", cell="", **far),),
    }


def _synthetic_claims() -> list[co.Claim]:
    far = ("chr1", 900_000, 900_400)
    return [
        claim(co.TARGET, "G1", cell="K562"),
        claim(co.ACTIVITY, "activates_target", gene="G1", cell="K562"),
        claim(co.CONTEXT, "K562", gene="G1", cell="K562"),
        claim(co.TARGET, "G2", cell="K562"),
        claim(co.ACTIVITY, "activates_target", gene="G2", cell="K562"),
        claim(co.ROLE, "enhancer_like"),
        claim(co.ORIGIN, "unique"),
        co.Claim("EH2", *far, co.TARGET, "G3", "G3", "K562"),
        co.Claim("EH2", *far, co.ACTIVITY, "activates_target", "G3", "K562"),
        co.Claim("EH2", *far, co.ROLE, "enhancer_like"),
        co.Claim("EH3", "chr2", 5, 500, co.TARGET, "G9", "G9", "placenta"),
    ]


def _report() -> co.Report:
    return co.judge(_synthetic_claims(), labels(), units=_synthetic_units(), references={})


def test_the_three_quantities_are_reported_apart_and_never_summed():
    r = _report()
    q = r.quantities()
    assert set(q) == {"target_accuracy", "role_accuracy", "coverage"}
    assert set(q["target_accuracy"]) == {co.TARGET}
    assert set(q["role_accuracy"]) == {co.ROLE, co.ACTIVITY}  # per axis, never pooled
    assert set(q["coverage"]) == set(co.AXES)
    shares = [s for block in q.values() for s in block.values()]
    with pytest.raises(TypeError):
        sum(shares)
    with pytest.raises(co.NotSummableError):
        q["target_accuracy"][co.TARGET] + q["coverage"][co.TARGET]
    with pytest.raises(co.NotSummableError):
        q["role_accuracy"][co.ROLE] + q["role_accuracy"][co.ACTIVITY]
    t = q["target_accuracy"][co.TARGET]
    assert (t.numerator, t.denominator) == (2, 3)  # G1 and G3 established, G2 refuted in K562
    cov = q["coverage"][co.TARGET]
    assert (cov.numerator, cov.denominator) == (3, 4)  # EH3 has no observation
    assert q["coverage"][co.ROLE].numerator == 0 and q["coverage"][co.ROLE].denominator == 2
    assert q["coverage"][co.ORIGIN].denominator == 1
    act = q["role_accuracy"][co.ACTIVITY]
    assert (act.numerator, act.denominator) == (1, 2)  # G1 fell (correct), G3 rose (incorrect)
    assert r.axes[co.ACTIVITY].not_judged_by_reason[co.CONDITION_UNMET] == 1  # G2 never responded


def _keys(x, out=None):
    out = set() if out is None else out
    if isinstance(x, dict):
        for k, v in x.items():
            out.add(str(k).lower())
            _keys(v, out)
    elif isinstance(x, list):
        for v in x:
            _keys(v, out)
    return out


def test_no_combined_score_anywhere_in_the_report():
    d = _report().to_dict()
    banned = {"score", "overall", "total", "combined", "mean", "average", "f1", "weighted", "summary_score"}
    assert not (_keys(d) & banned)
    assert d["no_sum"] == co.NO_SUM


def test_the_report_survives_serialisation():
    r = _report()
    d = json.loads(json.dumps(r.to_dict()))
    back = co.Report.from_dict(d)
    assert back.to_dict() == r.to_dict()
    for q, block in r.quantities().items():
        for a, s in block.items():
            assert co.Share.from_dict(d["quantities"][q][a]) == s


def test_heldout_file_pairs_judge_and_are_counted():
    r = _report()
    assert r.axes[co.TARGET].decided_with_heldout_file_pair == 1  # G3's increase is in the held-out split


# --- holdout's split discipline -----------------------------------------------------------------
def test_a_labelling_that_read_a_source_is_never_judged_by_it():
    units = _synthetic_units()
    with pytest.raises(ho.LeakError):
        co.judge(_synthetic_claims(), labels({"crispri:A"}), units=units, references={})
    with pytest.raises(ho.LeakError):
        co.judge(_synthetic_claims(), labels({ho.CRISPRI_HELDOUT_FILE}), units=units, references={})
    # built without S: judged by S alone; asking for every source raises
    rest_like = labels({"lentimpra:K562", "vista"}, built_without="crispri:A")
    r = co.judge(_synthetic_claims(), rest_like, units=units, sources=["crispri:A"], references={})
    assert r.sources == ["crispri:A"] and r.built_without == "crispri:A"
    with pytest.raises(ho.LeakError):
        co.judge(_synthetic_claims(), rest_like, units=units, references={})


def test_every_source_goes_through_holdouts_check(monkeypatch):
    seen = []
    real = ho.check_provenance

    def spy(lab, src, refs=None):
        seen.append(src)
        return real(lab, src, refs)

    monkeypatch.setattr(ho, "check_provenance", spy)
    co.judge(_synthetic_claims(), labels(), units=_synthetic_units(), references={})
    assert sorted(seen) == sorted(_synthetic_units())


def test_an_extra_observation_the_labelling_read_is_refused():
    extra = [
        co.Observation(
            "registry_biochemical",
            "encode_ccre_registry",
            value="enhancer_like",
            chrom=E[0],
            start=E[1],
            end=E[2],
        )
    ]
    with pytest.raises(ho.LeakError):
        co.judge(
            _synthetic_claims(), co.unchanged_labels(), units=_synthetic_units(), extra=extra, references={}
        )
    r = co.judge(_synthetic_claims(), labels(), units=_synthetic_units(), extra=extra, references={})
    assert r.quantities()["coverage"][co.ROLE].numerator == 1  # the registry judges EH1's role
    with pytest.raises(ValueError):
        co.judge([], labels(), units=_synthetic_units(), extra=[obs("crispri_decrease")], references={})


def test_the_unchanged_labels_read_no_holdout_source():
    lab = co.unchanged_labels()
    for s in ("crispri:Gasperini2019", "lentimpra:K562", "vista", "satmut", "gtex"):
        ho.check_provenance(lab, s, {})
    with pytest.raises(NotImplementedError):
        lab.predict(unit("vista", "positive"), "positive")


# --- the compiled labels as claims ---------------------------------------------------------------
PROGRAM = """\
region U_chr21_0 {
  locus: chr21:0-5000
  origin: assembly_gap
  molecular_role: unknown
  activity: unknown
  target_relation: unassigned
  evidence_status: curated_annotation
}
region U_chr21_5000 {
  locus: chr21:5000-9000
  origin: partly_repeat_derived/SINE, segmental_duplication
  molecular_role: promoter_like|unknown
  activity: unknown
  target_relation: unassigned
  evidence_status: curated_annotation, selection_not_detected
}
element EH1 {
  class: enhancer
  locus: chr21:10000-10400
  origin: unique
  molecular_role: enhancer_like, silencer|insulator_like|competing_promoter|unknown
  activity: represses_target
  target_relation: predicted_deletion_target
  evidence_status: registry_biochemical, predicted_model, selection_not_measured
}
{rule}
element EH1_measured {
  locus: chr21:10000-10400
  origin: unique
  molecular_role: enhancer_like
  activity: activates_target
  target_relation: measured_perturbation_target
  evidence_status: measured
}
rule EH1_measured activates G1 { strength: 0.5; when: cell_type = K562; evidence: experimental "CRISPRi" }
"""


RULE = 'rule EH1 inhibits G1 { strength: 0.2; when: cell_type = K562; evidence: predicted "AlphaGenome" }'


def test_compiled_claims(tmp_path):
    (tmp_path / "noncoding_chr21.bio").write_text(PROGRAM.replace("{rule}", RULE))
    nc: dict[str, Counter] = defaultdict(Counter)
    got = list(co.compiled_claims(tmp_path, not_claims=nc))
    by = Counter((c.element, c.axis, c.value) for c in got)
    assert by == Counter(
        {
            ("U_chr21_0", "origin", "assembly_gap"): 1,
            ("U_chr21_5000", "origin", "partly_repeat_derived/SINE"): 1,
            ("U_chr21_5000", "origin", "segmental_duplication"): 1,
            ("EH1", "origin", "unique"): 1,
            ("EH1", "molecular_role", "enhancer_like"): 1,
            ("EH1", "target", "G1"): 1,
            ("EH1", "activity", "represses_target"): 1,
            ("EH1", "context", "K562"): 1,
        }
    )
    assert not any(c.element.endswith("_measured") for c in got)  # the twins are evidence, not labels
    assert nc["molecular_role"] == Counter({"unknown": 1, "unchosen_alternatives": 2})
    assert nc["activity"]["region_unknown"] == 2 and nc["target"]["region_unassigned"] == 2
    t = next(c for c in got if c.axis == "target")
    assert (t.chrom, t.start, t.end, t.cell) == ("chr21", 10000, 10400, "K562")


def test_registration_carries_the_table():
    reg = co.registration()
    assert reg["registered"] == co.REGISTERED
    assert [r["observation"] for r in reg["table"]] == list(co.OBSERVATIONS)
    json.dumps(reg)
