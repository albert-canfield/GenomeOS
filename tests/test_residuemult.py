# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-residuemult: the lookup's definitions, the stated bound, and the guards on the record.

Every test here is predictor-side. None opens an RNA value, none fits anything, and none asserts a count
of real genes except through the committed record itself.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from genomeos.attribution import multrna as mr
from genomeos.attribution import residue as rs
from genomeos.attribution import residuemult as rm
from genomeos.attribution import saturate as sa
from genomeos.attribution import scalefit as sf

MODULE = Path("genomeos/attribution/residuemult.py")
SCRIPT = Path("scripts/residuemult_lookup.py")
RESULT_FILE = Path("data/results/residuemult_lookup.json")


# ==== nothing is declared twice and nothing is fitted ====================================================


def test_the_scale_is_scalefit_s_constant_and_is_not_re_declared_in_this_module() -> None:
    assert rm.SCALE == sf.SCALE
    assert rm.SCALE_ROUNDED == sf.SCALE_ROUNDED
    assert rm.SCALE == rs.SCALE
    body = MODULE.read_text()
    assert "0.2145911381525326" not in body, "the scale must be inherited, never re-typed"
    assert rm.CONVENTIONS["scale_is_refitted"] is False


def test_the_informativeness_threshold_is_7h_s_and_is_not_re_declared() -> None:
    assert rm.INFORMATIVE_LOG2 == mr.INFORMATIVE_LOG2
    assert "INFORMATIVE_LOG2 = multrna.INFORMATIVE_LOG2" in MODULE.read_text()


def test_no_module_or_script_of_this_lane_fits_optimises_or_re_estimates_anything() -> None:
    for path in (MODULE, SCRIPT):
        body = path.read_text()
        for word in ("curve_fit", "minimize", "least_squares", "fsolve", "brentq", "LinearRegression"):
            assert word not in body, f"{path} must not fit anything, found {word}"


def test_the_lane_declares_it_opens_no_rna_and_the_script_reads_only_the_symbol_column() -> None:
    body = SCRIPT.read_text()
    assert "NO RNA VALUE IS OPENED" in body
    # the only column of 7h's per-gene table the script names is `symbol`
    assert 'index("symbol")' in body
    assert '"D"' not in body and '"P"' not in body


# ==== the identity convention is stated, and it is stated for a reason ===================================


def test_the_match_field_is_the_symbol_and_the_record_says_why_the_identities_are_re_derived() -> None:
    assert rm.CONVENTIONS["match_field"] == "symbol"
    assert "`symbol`" in rm.IDENTITY
    assert "carries NEITHER" in rm.IDENTITY
    assert "RE-DERIVED" in rm.IDENTITY
    assert "combine2.gene_rules" in rm.RE_DERIVED
    assert "scalefit.scaled_rho" in rm.RE_DERIVED


def test_lane_residue_s_committed_per_gene_table_really_carries_no_gene_identity() -> None:
    """The fact IDENTITY rests on: if this ever changes, the re-derivation is no longer necessary."""
    import json

    import tracked_paths as tp

    tp.must_be_present(
        Path("data/results/residue_genes.json"),
        was="lane-residue's committed per-gene table, 5d2c8f5",
    )
    cols = json.loads(Path("data/results/residue_genes.json").read_text())["columns"]
    assert not {"gene", "symbol", "gene_id"} & set(cols), cols
    assert "chrom" in cols


# ==== the multiplicative form has no basal share, which is the whole point ===============================


def test_the_multiplicative_form_cannot_refuse_a_gene_for_basal_not_positive() -> None:
    """Why the residue is the ADDITIVE form's failure: the reason it was excluded for has no state here."""
    from genomeos.attribution import combine2

    assert "basal_not_positive" not in rm.saturate_refusals()
    # a gene whose drops sum far past its level: the mean refuses it, the multiplicative form does not
    rho = {f"e{i}": 0.5 for i in range(10)}
    mean = combine2.fit("mean", rho)
    assert not mean.ok and mean.reason == "basal_not_positive"
    assert sa.fit(rho).ok


def test_every_declared_refusal_is_one_saturate_can_actually_give() -> None:
    body = Path("genomeos/attribution/saturate.py").read_text()
    for reason in rm.saturate_refusals():
        assert f'"{reason}"' in body, reason


def test_the_scale_relaxes_the_inhibitor_bound_and_never_tightens_it() -> None:
    """WHY_THE_SCALE_CAN_ONLY_HELP's runnability half, as arithmetic on one gene at a time."""
    # RHO_j = 3 is above the unscaled bound of 2 and below the scaled bound of 1 + 1/SCALE
    rho = {"a": 0.4, "j": 3.0}
    assert not sa.fit(rho).ok
    assert sa.fit(rho).reason == "inhibitor_strength_above_1"
    assert sa.fit(sf.scaled_rho(rho, rm.SCALE)).ok
    assert pytest.approx(5.6601, abs=1e-4) == 1.0 + 1.0 / rm.SCALE


def test_scaling_cannot_take_runnability_away_from_any_gene_the_form_already_ran() -> None:
    for rho in (
        {"a": 0.4, "b": 0.2, "c": 0.9},
        {"a": 1e-8, "b": 0.5},
        {"a": 0.3, "j": 1.9},
        {"a": 0.3, "j": 1.999999},
    ):
        if sa.fit(rho).ok:
            assert sa.fit(sf.scaled_rho(rho, rm.SCALE)).ok, rho


def test_the_record_states_the_informativeness_caveat_rather_than_claiming_a_theorem() -> None:
    assert "ONE CAVEAT" in rm.WHY_THE_SCALE_CAN_ONLY_HELP_THE_MULTIPLICATIVE_FORM
    assert "near-certainty, not a theorem" in rm.WHY_THE_SCALE_CAN_ONLY_HELP_THE_MULTIPLICATIVE_FORM


# ==== the bound checks catch what they say they catch =====================================================


class _FakeGeneRules:
    """A gene with a reader state per cell, enough for the multiplicative closed form's range."""

    def __init__(self, rho: dict[str, float], open_by_cell: list[dict[str, bool]]) -> None:
        self.rho = rho
        self.symbol = "FAKE1"
        self.chrom = "chrZ"
        self.elements = list(rho)
        self._open = open_by_cell

    def open_in(self, j: int) -> dict[str, bool]:
        return self._open[j]


def _fake(rho: dict[str, float], vary: bool) -> _FakeGeneRules:
    n = len(mr.CELLS)
    if vary:
        opens = [{el: (j % 2 == 0) for el in rho} for j in range(n)]
    else:
        opens = [dict.fromkeys(rho, True) for _ in range(n)]
    return _FakeGeneRules(rho, opens)


def test_informative_is_runnable_and_a_level_that_differs_between_reader_cells() -> None:
    rho = {"a": 0.4, "b": 0.7}
    inf, rng = rm.informative_under(sa.fit(rho), _fake(rho, vary=True))
    assert inf and rng > rm.INFORMATIVE_LOG2
    inf, rng = rm.informative_under(sa.fit(rho), _fake(rho, vary=False))
    assert not inf and rng == 0.0


def test_a_gene_the_form_cannot_run_is_never_informative_and_its_range_is_not_a_number() -> None:
    rho = {"a": 0.4, "j": 3.0}
    inf, rng = rm.informative_under(sa.fit(rho), _fake(rho, vary=True))
    assert not inf and math.isnan(rng)


def test_PLANTED_an_undeclared_refusal_is_reported_as_a_bound_violation() -> None:
    row = rm.ResidueGene(
        symbol="PLANT1",
        chrom="chrZ",
        n_rules=3,
        n_inhibiting=0,
        in_committed_table=False,
        runnable_unscaled=False,
        informative_unscaled=False,
        runnable_scaled=True,
        informative_scaled=True,
        reason_unscaled="basal_not_positive",  # a reason the multiplicative form has no state for
        reason_scaled="",
        range_log2_scaled=1.0,
    )
    bad = rm.bound_violations(row)
    assert any("not a declared one" in b for b in bad), bad


def test_PLANTED_a_gene_runnable_unscaled_but_not_at_the_scale_is_reported() -> None:
    row = rm.ResidueGene(
        symbol="PLANT2",
        chrom="chrZ",
        n_rules=3,
        n_inhibiting=1,
        in_committed_table=True,
        runnable_unscaled=True,
        informative_unscaled=True,
        runnable_scaled=False,
        informative_scaled=False,
        reason_unscaled="",
        reason_scaled="inhibitor_strength_above_1",
        range_log2_scaled=float("nan"),
    )
    bad = rm.bound_violations(row)
    assert any("runnable unscaled but not at the scale" in b for b in bad), bad
    assert any("informative unscaled but not at the scale" in b for b in bad), bad


def test_a_clean_row_reports_no_violation_so_the_check_is_not_vacuous() -> None:
    row = rm.ResidueGene(
        symbol="CLEAN1",
        chrom="chrZ",
        n_rules=4,
        n_inhibiting=0,
        in_committed_table=True,
        runnable_unscaled=True,
        informative_unscaled=True,
        runnable_scaled=True,
        informative_scaled=True,
        reason_unscaled="",
        reason_scaled="",
        range_log2_scaled=0.5,
    )
    assert rm.bound_violations(row) == []
    assert row.fitted


def test_the_union_a_row_reports_is_count_1_or_count_2_and_nothing_else() -> None:
    def row(in_table: bool, inf_scaled: bool) -> rm.ResidueGene:
        return rm.ResidueGene(
            symbol="U",
            chrom="chrZ",
            n_rules=2,
            n_inhibiting=0,
            in_committed_table=in_table,
            runnable_unscaled=True,
            informative_unscaled=True,
            runnable_scaled=inf_scaled,
            informative_scaled=inf_scaled,
            reason_unscaled="",
            reason_scaled="",
            range_log2_scaled=1.0,
        )

    assert row(True, True).fitted
    assert row(True, False).fitted
    assert row(False, True).fitted
    assert not row(False, False).fitted


# ==== the implication rule, its threshold and its wording ================================================


def test_both_branches_of_the_implication_rule_are_reachable_at_the_declared_threshold() -> None:
    assert rm.MOST == 0.5
    assert rm.implication(1854, 3708) == "NOT_ANSWERED"  # exactly a half is NOT "most"
    assert rm.implication(1855, 3708) == "ALREADY_ANSWERED"  # one gene more is
    assert rm.implication(3708, 3708) == "ALREADY_ANSWERED"
    assert rm.implication(0, 3708) == "NOT_ANSWERED"


def test_exactly_a_half_is_not_most_and_one_gene_more_is() -> None:
    assert rm.implication(500, 1000) == "NOT_ANSWERED"
    assert rm.implication(501, 1000) == "ALREADY_ANSWERED"


def test_an_empty_residue_is_refused_rather_than_decided() -> None:
    with pytest.raises(ValueError):
        rm.implication(0, 0)


def test_the_already_answered_reading_says_no_new_form_is_needed_in_those_words() -> None:
    text = rm.READINGS["ALREADY_ANSWERED"]
    assert "ADDITIVE FORM'S STRUCTURAL FAILURE" in text
    assert "ALREADY ANSWERED BY THE FORM M2h VALIDATED" in text
    assert "NO NEW FORM IS NEEDED" in text


def test_the_not_answered_reading_proposes_nothing() -> None:
    text = rm.READINGS["NOT_ANSWERED"].lower()
    assert "proposes nothing" in text
    for word in ("cap", "saturat", "recommend", "next form"):
        assert word not in text, word


def test_every_branch_of_the_rule_has_a_reading_and_no_reading_lacks_a_branch() -> None:
    assert set(rm.READINGS) == {"ALREADY_ANSWERED", "NOT_ANSWERED"}
    for branch in rm.READINGS:
        text = rm.reading(branch, union=1, residue=2, share="50.0%", committed=1, scaled=1, complement=1)
        assert "{" not in text, text


# ==== what the lane refuses to claim =====================================================================


def test_the_module_carries_the_theorem_that_makes_a_saturating_census_pointless() -> None:
    assert "WHY_NOT_SATURATING" in MODULE.read_text()
    assert "concave" in rm.__doc__ and "CONVEX" in rm.__doc__
    assert "redundancy" in rm.__doc__


def test_the_lane_carries_multrna_s_worse_outcome_forward_word_for_word() -> None:
    joined = " ".join(rm.LIMITATIONS)
    assert "outcome on measured RNA was WORSE" in joined
    assert "lane-recovered" in joined
    assert "a fit accepts a gene" in joined


def test_the_module_discloses_its_prior_exposure_including_the_bound_it_already_knew() -> None:
    assert "11,396" in rm.PRIOR_EXPOSURE
    assert "at most 1,104" in rm.PRIOR_EXPOSURE
    assert "surprise" in rm.PRIOR_EXPOSURE


def test_the_module_claims_nothing_about_expression_and_registers_no_test() -> None:
    doc = rm.__doc__
    assert "reads NO RNA VALUE" in doc
    assert "registers NO TEST" in doc
    assert "no claim about expression" in doc


def test_the_stop_list_check_is_named_in_the_module() -> None:
    doc = rm.__doc__
    assert "STOP-LIST CHECK" in doc
    for item in ("AlphaGenome", "deletion", "biosample", "area J", "node-model", "BioForge"):
        assert item in doc, item


def test_every_bound_constant_is_a_real_attribute_so_the_record_cannot_drift() -> None:
    for field, const in rm.BOUND.items():
        assert hasattr(rm, const), (field, const)


# ==== the committed record: the counts are RECOMPUTED from the table, never read as written ==============


def _record() -> dict:
    import json

    import tracked_paths as tp

    tp.must_be_present(
        RESULT_FILE,
        was="the one lookup run of this lane, written at the committed definitions",
    )
    return json.loads(RESULT_FILE.read_text())


def test_every_published_count_recomputes_from_the_per_gene_table() -> None:
    """A committed number must be recomputable from the rows beside it (the lesson of 2026-10-12)."""
    d = _record()
    t = d["per_gene_table"]
    i = {c: t["columns"].index(c) for c in t["columns"]}
    rows = t["rows"]
    c = d["counts"]
    assert len(rows) == d["residue_genes"]
    assert sum(1 for r in rows if r[i["in_committed_table"]]) == c["count_1_in_multrna_s_committed_table"]
    assert (
        sum(1 for r in rows if r[i["informative_unscaled"]])
        == c["count_1b_informative_multiplicatively_predictor_side"]
    )
    assert (
        sum(1 for r in rows if r[i["informative_scaled"]])
        == c["count_2_informative_multiplicatively_at_the_scale"]
    )
    assert (
        sum(1 for r in rows if r[i["in_committed_table"]] or r[i["informative_scaled"]])
        == c["union_of_count_1_and_count_2"]
    )
    assert (
        sum(1 for r in rows if not (r[i["in_committed_table"]] or r[i["informative_scaled"]]))
        == c["informative_under_neither"]
    )


def test_the_two_counts_and_the_complement_close_on_the_population() -> None:
    d = _record()
    c = d["counts"]
    assert c["union_of_count_1_and_count_2"] + c["informative_under_neither"] == d["residue_genes"]
    assert c["count_1_in_multrna_s_committed_table"] <= c["union_of_count_1_and_count_2"]
    assert c["count_2_informative_multiplicatively_at_the_scale"] <= c["union_of_count_1_and_count_2"]


def test_the_branch_and_the_implication_recompute_by_the_committed_rule() -> None:
    d = _record()
    branch = rm.implication(d["counts"]["union_of_count_1_and_count_2"], d["residue_genes"])
    assert branch == d["branch"]
    assert d["most"] == rm.MOST
    head = rm.READINGS[branch].split("{")[0]
    assert d["implication"].startswith(head)


def test_the_re_derived_population_reconciles_with_lane_residue_s_committed_counts() -> None:
    """The re-derivation is the lane's own run: it must land on lane-residue's 3,708, or say it did not."""
    d = _record()
    rec = d["reconciliation"]
    assert rec["agrees"] is (rec["differences"] == {})
    assert rec["re_derived"]["the_residue"] == d["residue_genes"]
    assert rec["agrees"], rec["differences"]
    assert rec["committed"]["excluded_basal_not_positive"] == (
        rec["committed"]["recovered_at_the_scale"] + rec["committed"]["the_residue"]
    )


def test_the_stated_bound_held_on_every_gene_and_no_violation_was_absorbed() -> None:
    d = _record()
    assert d["bound_violations"] == []
    c = d["counts"]
    # runnability cannot be lost to the scale: the exact half of the stated bound
    assert c["runnable_at_the_scale"] >= c["runnable_unscaled"]
    assert (
        c["count_2_informative_multiplicatively_at_the_scale"]
        >= c["count_1b_informative_multiplicatively_predictor_side"]
    )
    assert d["identity_checks"]["max_abs_q_after_minus_scale_times_q_before"] < 1e-9


def test_no_residue_gene_is_refused_multiplicatively_for_the_reason_the_mean_refused_it() -> None:
    """The structural point, checked on the real population and not only on a constructed gene."""
    d = _record()
    for where in ("unscaled", "at_the_scale"):
        reasons = set(d["reasons_not_runnable"][where])
        assert "basal_not_positive" not in reasons, (where, reasons)
        assert reasons <= rm.saturate_refusals(), (where, reasons)


def test_the_complement_is_reported_by_reason_so_it_cannot_read_as_genes_the_form_fits() -> None:
    d = _record()
    by_reason = d["informative_under_neither_by_reason"]
    assert sum(by_reason.values()) == d["counts"]["informative_under_neither"]
    assert set(by_reason) <= rm.saturate_refusals() | {"runnable_not_informative"}


def test_the_record_carries_the_match_field_the_limitations_and_the_refusal_to_claim_expression() -> None:
    d = _record()
    assert d["match_field"] == rm.CONVENTIONS["match_field"] == "symbol"
    assert d["rna_read"] == "none"
    assert d["alphagenome_requests"] == 0
    assert "NO RNA VALUE WAS OPENED IN THIS LANE" in d["what_this_is_not"]
    assert "WORSE" in d["what_this_is_not"]
    assert d["limitations"] == list(rm.LIMITATIONS)
    assert d["prior_exposure"] == rm.PRIOR_EXPOSURE
    assert d["implication_rule"] == rm.IMPLICATION_RULE


def test_the_record_reports_whether_a_residue_symbol_names_more_than_one_gene_in_this_run() -> None:
    """LIMITATIONS names symbol ambiguity as the match's weakness; the record must carry the count."""
    d = _record()
    assert isinstance(d["identity_checks"]["residue_symbols_not_unique_in_this_run"], list)
    t = d["per_gene_table"]
    i = t["columns"].index("symbol")
    syms = [r[i] for r in t["rows"]]
    assert len(set(syms)) + len(d["identity_checks"]["residue_symbols_not_unique_in_this_run"]) >= len(syms)
