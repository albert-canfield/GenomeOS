# SPDX-License-Identifier: AGPL-3.0-or-later
"""The dependency a deep deletion creates, wired into the therapeutics pipeline's evidence.

Registered 2026-10-10 in docs/THERAPEUTICS.md ("The dependency a deep deletion
creates"). These tests hold the registration's conditions: the trigger is a deep
deletion and nothing else, the partners are the sweep's and no others, the lineage
verdict is unknown wherever it cannot be measured, the evidence is never the
patient's, and no target candidate moves.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.cancer.alterations import GeneAlteration, grade
from genomeos.therapeutics import dependency as dep
from genomeos.therapeutics import machine_report, text_report
from genomeos.therapeutics import synthetic_lethality as sl
from genomeos.therapeutics.design import dataset
from genomeos.therapeutics.model import TherapeuticTargetCandidate
from genomeos.therapeutics.pipeline import PatientProfile, alteration_origin, analyse, limitations_for
from tests.test_therapeutics import LOW, StubProviders, hpa_row

SWEEP = Path("data/results/synthetic_lethality.json")
FORBIDDEN = ("effective", "will respond", "dose", "dosing", "treat", "cure", "protocol")


def deleted(*genes: str) -> list[GeneAlteration]:
    return grade([GeneAlteration(g, "deep_deletion", gistic=-2) for g in genes], knowledge={})


def view(
    deleted_lines: int = 8,
    kept_lines: int = 30,
    deleted_effect: float = -1.0,
    lineage: str = "CNS/Brain",
    mtap_kept_null_lines: int = 0,
    spread: float | None = None,
) -> dep.DepMapView:
    """A small DepMap: lines that deleted CDKN2A (with MTAP) depend on WDR77 and PRMT5.

    `mtap_kept_null_lines` adds CDKN2A-deleted lines that kept MTAP and carry no
    dependency at all, interleaved with the kept lines: the same-event control.
    `spread` places the deleted lines' effects at `deleted_effect + spread * i` for
    both partners, so they can be interleaved with the kept lines' 0.000 to 0.029.
    """
    effect: dict[str, dict[str, float]] = {}
    copies: dict[str, dict[str, float]] = {}
    lineage_of: dict[str, str] = {}
    for i in range(deleted_lines):
        m = f"DEL{i}"
        effect[m] = (
            {"WDR77": deleted_effect - 0.01 * i, "PRMT5": deleted_effect - 0.02 * i}
            if spread is None
            else {"WDR77": deleted_effect + spread * i, "PRMT5": deleted_effect + spread * i}
        )
        copies[m] = {"CDKN2A": 0.0, "MTAP": 0.1, "CDKN2B": 0.0}
        lineage_of[m] = lineage
    for i in range(mtap_kept_null_lines):
        m = f"NUL{i}"
        effect[m] = {"WDR77": 0.0021 + 0.005 * i, "PRMT5": 0.0021 + 0.005 * i}
        copies[m] = {"CDKN2A": 0.0, "MTAP": 2.0, "CDKN2B": 0.0}
        lineage_of[m] = lineage
    for i in range(kept_lines):
        m = f"KEP{i}"
        effect[m] = {"WDR77": 0.001 * i, "PRMT5": 0.001 * i}
        copies[m] = {"CDKN2A": 2.0, "MTAP": 2.0, "CDKN2B": 2.0}
        lineage_of[m] = lineage
    # another lineage whose lines would flip the answer if they were read
    for i in range(10):
        m = f"OTH{i}"
        effect[m] = {"WDR77": 1.0, "PRMT5": 1.0}
        copies[m] = {"CDKN2A": 0.0, "MTAP": 0.0, "CDKN2B": 0.0}
        lineage_of[m] = "Lung"
    return dep.DepMapView(effect, copies, frozenset({"CDKN2A", "CDKN2B", "MTAP"}), lineage_of, True, "")


# --- the frozen partner table is the sweep's ---------------------------------------------


def test_the_registered_pairs_are_the_committed_sweep_rows():
    committed = json.loads(SWEEP.read_text())
    assert committed["hits"] == len(dep.REGISTERED_PAIRS) == 44
    for hit, pair in zip(committed["top_hits"], dep.REGISTERED_PAIRS, strict=False):
        a, b, n_lost, n_kept, m_lost, m_kept, delta, q = pair
        assert (hit["lost"], hit["dependency"]) == (a, b)
        assert (hit["lines_lost"], hit["lines_intact"]) == (n_lost, n_kept)
        assert (hit["median_effect_lost"], hit["median_effect_intact"], hit["cliffs_delta"]) == (
            m_lost,
            m_kept,
            delta,
        )
        assert abs(hit["q"] - q) <= 5e-6, "the committed q is rounded to five decimals"
    assert all(q <= sl.FDR for *_, q in dep.REGISTERED_PAIRS)
    assert {tuple(e["genes"]) for e in committed["loss_genes_that_are_the_same_event"]} == {
        (a, b) for a, b, _ in dep.SAME_EVENT
    }
    assert len(dep.SWEEP_TESTABLE) == committed["loss_genes_with_enough_lines"] == 171
    assert len(dep.SWEEP_ASKED) == 304 and not dep.SWEEP_TESTABLE & dep.SWEEP_ASKED_UNTESTABLE
    assert {a for a, *_ in dep.REGISTERED_PAIRS} <= dep.SWEEP_TESTABLE


@pytest.mark.needs_local_data(
    "data/knowledge/therapeutics/depmap", how="the synthetic-lethality sweep's cache"
)
def test_the_frozen_table_is_what_the_sweep_returns_from_its_cache():
    """The falsifier: a re-run of the sweep offline gives back exactly the frozen table."""
    from genomeos.therapeutics import depmap

    paralogues = sl.paralogue_table(list(sl.PANEL), net=False)
    pairs, _space = sl.candidate_pairs(paralogues)
    dependencies = sorted({b for _a, b in pairs})
    loss_genes = sorted({a for a, _b in pairs} | {c["lost"] for c in sl.PREREGISTERED if c["lost"] != "MSI"})
    effect = depmap.matrix("gene_effect", set(dependencies), net=False)["values"]
    damaging = depmap.matrix("damaging", set(loss_genes), net=False)["values"]
    copies = depmap.matrix("absolute_cn", set(loss_genes), net=False)["values"]
    sig = depmap.signatures(net=False)
    msi = {m: v for m, v in ((m, sl._float(r.get("MSIScore"))) for m, r in sig.items()) if v is not None}
    store = sl._RankStore(effect)
    statuses = {}
    for gene in loss_genes:
        status = sl.loss_status(gene, "loss", damaging, copies, msi)
        if sum(1 for v in status.values() if v) >= sl.MIN_LOST:
            statuses[gene] = status
    coverage = {g: store.canonical(frozenset(s)) for g, s in statuses.items()}
    rows = [sl.test_pair(a, b, "loss", statuses[a], store, coverage[a]) for a, b in pairs if a in statuses]
    tested = [r for r in rows if "p" in r]
    for r, q in zip(tested, sl.benjamini_hochberg([r["p"] for r in tested]), strict=True):
        r["q"] = q
    hits = sorted(
        (r for r in tested if r["q"] <= sl.FDR and r["cliffs_delta"] < 0),
        key=lambda r: (r["q"], r["cliffs_delta"]),
    )
    assert len(tested) == 19749
    assert [(r["lost"], r["dependency"], r["lines_lost"], r["cliffs_delta"]) for r in hits] == [
        (a, b, n, d) for a, b, n, _k, _ml, _mk, d, _q in dep.REGISTERED_PAIRS
    ]
    assert set(loss_genes) == dep.SWEEP_ASKED
    assert set(statuses) == dep.SWEEP_TESTABLE


# --- the trigger --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "event",
    [
        GeneAlteration("CDKN2A", "shallow_deletion", gistic=-1),
        GeneAlteration("CDKN2A", "amplification", copies=9.0),
        GeneAlteration("CDKN2A", "gain", gistic=1),
        GeneAlteration("ARID1A", "fusion", partner="X"),
    ],
)
def test_only_a_deep_deletion_creates_a_dependency(event):
    out = dep.assess(grade([event], knowledge={}), {"CDKN2A": 1.0, "MTAP": 2.0}, "CNS/Brain", view=view())
    assert out == dep.empty()


def test_a_shallow_deletion_beside_a_deep_one_adds_nothing_of_its_own():
    alterations = grade(
        [
            GeneAlteration("CDKN2A", "deep_deletion", gistic=-2),
            GeneAlteration("ARID1A", "shallow_deletion", gistic=-1),
        ],
        knowledge={},
    )
    out = dep.assess(alterations, {}, "CNS/Brain", view=view())
    assert [r["gene"] for r in out["dependencies"]] == ["CDKN2A"]
    assert {c["lost_gene"] for c in out["dependency_candidates"]} == {"CDKN2A"}


@pytest.mark.parametrize(
    "p,q,delta,verdict",
    [
        (0.04, 0.09, -0.3, "supported_in_lineage"),
        (0.04, 0.2, -0.3, "nominal_in_lineage"),
        (0.04, 0.2, 0.3, "not_supported_in_lineage"),
        (0.001, 0.001, 0.2, "not_supported_in_lineage"),
        (0.2, 0.3, -0.3, "not_supported_in_lineage"),
    ],
)
def test_the_verdict_rule_is_the_sweeps(p, q, delta, verdict):
    assert dep.lineage_verdict(p, q, delta) == verdict


def test_a_diploid_count_is_not_a_deletion():
    assert dep.assess([], {"CDKN2A": 2.0, "MTAP": 0.9}, "CNS/Brain", view=view()) == dep.empty()


@pytest.mark.parametrize(
    "gene,status",
    [("PTEN", "not_asked"), ("RB1", "not_asked"), ("PRMT5", "asked_untestable"), ("ATR", "asked_no_partner")],
)
def test_a_deletion_without_a_registered_partner_says_why_and_adds_nothing(gene, status):
    out = dep.assess(deleted(gene), {}, "CNS/Brain", view=view())
    (record,) = out["dependencies"]
    assert record["status"] == status and record["partners"] == []
    assert out["dependency_candidates"] == []


# --- the lineage verdicts -----------------------------------------------------------------


def test_unstated_lineage_is_unknown_and_reads_no_cell_lines():
    out = dep.assess(deleted("CDKN2A"), {}, "", ["CDKN2A"])
    assert [c["gene"] for c in out["dependency_candidates"]] == ["WDR77", "PRMT5"]
    for c in out["dependency_candidates"]:
        assert c["lineage"]["verdict"] == "unknown"
        assert "not stated" in c["lineage"]["reason"]
        assert "p" not in c["lineage"]
    assert out["dependency_lineage"]["read"] is None


def test_an_unrecognised_lineage_is_unknown():
    out = dep.assess(deleted("CDKN2A"), {}, "Atlantis", view=view())
    assert {c["lineage"]["verdict"] for c in out["dependency_candidates"]} == {"unknown"}
    assert "not a DepMap OncotreeLineage" in out["dependency_candidates"][0]["lineage"]["reason"]


def test_the_lineage_test_reads_only_the_stated_lineage_and_is_case_blind():
    out = dep.assess(deleted("CDKN2A"), {}, "cns/brain", view=view())
    for c in out["dependency_candidates"]:
        t = c["lineage"]
        assert t["name"] == "CNS/Brain"
        assert t["verdict"] == "supported_in_lineage", t
        assert (t["lines_lost"], t["lines_intact"]) == (8, 30), "the Lung lines would have flipped it"
        assert t["tests_in_this_analysis"] == 2


@pytest.mark.parametrize("deleted_lines,kept_lines", [(4, 30), (8, 19)])
def test_too_few_lines_is_unknown_with_the_counts(deleted_lines, kept_lines):
    out = dep.assess(deleted("CDKN2A"), {}, "CNS/Brain", view=view(deleted_lines, kept_lines))
    for c in out["dependency_candidates"]:
        t = c["lineage"]
        assert t["verdict"] == "unknown"
        assert (t["lines_lost"], t["lines_intact"]) == (deleted_lines, kept_lines)
        assert f"at least {sl.MIN_LOST} deleted and {sl.MIN_INTACT} kept" in t["reason"]


def test_no_shift_is_not_supported():
    out = dep.assess(deleted("CDKN2A"), {}, "CNS/Brain", view=view(deleted_effect=0.0005, spread=0.004))
    assert {c["lineage"]["verdict"] for c in out["dependency_candidates"]} == {"not_supported_in_lineage"}


def test_a_partner_without_a_gene_effect_and_a_masked_lost_gene_are_unknown():
    v = view()
    for row in v.effect.values():
        row.pop("PRMT5")
    out = dep.assess(deleted("CDKN2A"), {}, "CNS/Brain", view=v)
    prmt5 = next(c for c in out["dependency_candidates"] if c["gene"] == "PRMT5")
    assert prmt5["lineage"]["verdict"] == "unknown" and "no gene effect" in prmt5["lineage"]["reason"]
    masked = dep.DepMapView(v.effect, v.copies, frozenset({"MTAP"}), v.lineage_of, True, "")
    out = dep.assess(deleted("CDKN2A"), {}, "CNS/Brain", view=masked)
    assert all("no absolute copy number" in c["lineage"]["reason"] for c in out["dependency_candidates"])


def test_an_absent_cache_is_unknown_and_the_sweep_rows_still_stand():
    absent = dep.DepMapView(reason="the DepMap cache is not available on this machine (test)")
    out = dep.assess(deleted("CDKN2A"), {}, "CNS/Brain", view=absent)
    for c in out["dependency_candidates"]:
        assert c["lineage"]["verdict"] == "unknown" and "not available" in c["lineage"]["reason"]
        assert c["sweep"]["lines_lost"] == 427


def test_a_partner_this_tumour_also_deleted_is_contradicted_and_listed_last():
    out = dep.assess(deleted("CDKN2A", "WDR77"), {}, "CNS/Brain", view=view())
    verdicts = [(c["gene"], c["lineage"]["verdict"]) for c in out["dependency_candidates"]]
    assert verdicts[-1] == ("WDR77", "contradicted_by_this_tumour")
    assert out["dependency_candidates"][-1]["partner_state_in_this_tumour"] == "deleted"


def test_the_same_event_gene_restricts_the_deleted_lines_to_ones_like_this_tumour():
    """The control: CDKN2A-deleted lines that kept MTAP carry no dependency here."""
    v = view(deleted_lines=8, mtap_kept_null_lines=8)
    like_mtap_kept = dep.assess(deleted("CDKN2A"), {"MTAP": 2.0}, "CNS/Brain", view=v)
    like_mtap_deleted = dep.assess(deleted("CDKN2A", "MTAP"), {}, "CNS/Brain", view=v)
    mixed = dep.assess(deleted("CDKN2A"), {}, "CNS/Brain", view=v)
    kept = {
        c["gene"]: c["lineage"] for c in like_mtap_kept["dependency_candidates"] if c["lost_gene"] == "CDKN2A"
    }
    assert {t["verdict"] for t in kept.values()} == {"not_supported_in_lineage"}
    assert all(t["lines_lost"] == 8 for t in kept.values())
    assert kept["WDR77"]["restricted_to_lines_matching_this_tumour_on"] == [{"gene": "MTAP", "state": "kept"}]
    co = [c for c in like_mtap_deleted["dependency_candidates"] if c["lost_gene"] == "CDKN2A"]
    assert {c["lineage"]["verdict"] for c in co} == {"supported_in_lineage"}
    assert all(c["lineage"]["lines_lost"] == 8 for c in co)
    for c in mixed["dependency_candidates"]:
        assert c["lineage"]["lines_lost"] == 16
        assert c["lineage"]["restricted_to_lines_matching_this_tumour_on"] == []
        assert "mixes lines that lost it and lines that kept it" in c["limitations"][0]


# --- the evidence and its wording ---------------------------------------------------------


def test_dependency_evidence_is_never_the_patients_and_never_above_seven_tenths():
    out = dep.assess(deleted("CDKN2A", "MTAP"), {}, "CNS/Brain", view=view())
    assert out["dependency_candidates"]
    for c in out["dependency_candidates"]:
        assert len(c["evidence"]) == (2 if "p" in c["lineage"] else 1)
        for e in c["evidence"]:
            assert e["patient_derived"] is False
            assert e["level"] in ("computational", "in_vitro")
            assert e["confidence"] <= 0.7
    for r in out["dependencies"]:
        (e,) = r["evidence"]
        assert e["patient_derived"] is True and r["gene"] in e["claim"]


def test_every_reading_ends_with_the_hypothesis_sentence_and_claims_no_therapy():
    v = view()
    outs = [
        dep.assess(deleted("CDKN2A", "WDR77"), {}, "CNS/Brain", view=v),
        dep.assess(deleted("CDKN2A"), {}, "", ["CDKN2A"]),
        dep.assess(deleted("CDKN2A"), {}, "CNS/Brain", view=view(deleted_effect=0.0005, spread=0.004)),
        dep.assess(deleted("CDKN2A"), {}, "CNS/Brain", view=view(4, 30)),
    ]
    seen = set()
    for out in outs:
        for c in out["dependency_candidates"]:
            seen.add(c["lineage"]["verdict"])
            assert c["reading"].endswith(dep.CLOSING)
            low = (c["reading"] + " " + c["route"] + " " + " ".join(c["limitations"])).lower()
            assert not [w for w in FORBIDDEN if w in low], c["reading"]
    assert seen == {
        "supported_in_lineage",
        "unknown",
        "not_supported_in_lineage",
        "contradicted_by_this_tumour",
    }


def test_the_order_is_verdict_then_the_sweeps_q():
    out = dep.assess(deleted("MTAP", "CDKN2A"), {}, "CNS/Brain", view=view())
    keys = [
        (dep.VERDICT_ORDER.index(c["lineage"]["verdict"]), c["sweep"]["q"])
        for c in out["dependency_candidates"]
    ]
    assert keys == sorted(keys)
    assert {c["gene"] for c in out["dependency_candidates"]} == {"WDR77", "PRMT5", "MAT2A"}


# --- the pipeline: candidates added, targets untouched ------------------------------------


def _profile() -> PatientProfile:
    alterations = grade(
        [
            GeneAlteration("ERBB2", "amplification", copies=12.0),
            GeneAlteration("CDKN2A", "deep_deletion", copies=0.0),
        ],
        knowledge={},
    )
    return PatientProfile(
        "t",
        "t.vcf",
        copy_number={"ERBB2": 12.0, "CDKN2A": 0.0},
        alterations=alterations,
        alterations_path="t.cnv",
    )


def _targets(a: dict) -> list[tuple]:
    return [
        (c.gene, c.scores.overall, c.evidence_tier, c.mechanism_reach, c.target_class, c.best_mechanism)
        for c in a["candidates"]
    ]


def test_the_pipeline_adds_dependency_candidates_and_moves_no_target(monkeypatch):
    providers = StubProviders(expression={"ERBB2": hpa_row("ERBB2", LOW)})
    on = analyse([], _profile(), providers)
    monkeypatch.setattr(dep, "assess", lambda *a, **k: dep.empty())
    off = analyse([], _profile(), providers)
    assert _targets(on) == _targets(off), "a target candidate moved"
    assert off["dependency_candidates"] == []
    added = [c["gene"] for c in on["dependency_candidates"]]
    assert added == ["WDR77", "PRMT5"]
    assert not set(added) & {c.gene for c in on["candidates"]}, "a dependency entered the target ranking"
    assert on["dependencies"][0]["gene"] == "CDKN2A"
    report = machine_report(on)
    assert [c["gene"] for c in report["dependency_candidates"]] == added
    assert not set(added) & {c["gene"] for c in report["therapeutic_candidates"]}
    assert "DEPENDENCIES A DEEP DELETION CREATES" in text_report(on)
    assert "DEPENDENCIES A DEEP DELETION CREATES" not in text_report(off)
    design = json.dumps(dataset(on, providers))
    assert "WDR77" not in design and "PRMT5" not in design


def test_the_deleted_genes_limitation_points_at_the_list():
    c = TherapeuticTargetCandidate(gene="CDKN2A")
    c.origins = [alteration_origin(GeneAlteration("CDKN2A", "deep_deletion", gistic=-2), "t")]
    (line,) = [x for x in limitations_for(c) if "homozygous deletion" in x]
    assert "dependency_candidates" in line and "does not model" not in line
