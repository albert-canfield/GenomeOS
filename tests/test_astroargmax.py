# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item (h): what the astrocyte-argmax REGISTRATION must be, tested before any answer is bought.

The checks that matter here are about ORDER and about what is ABSENT:

* the bands exist as numbers and are RATIOS, so no value of the control can have fixed them;
* the control is a share of COMPILED RULES and never of track names, joined by CURIE;
* the panel share is labelled a panel share everywhere it appears;
* no field a run would fill exists in the payload, not even empty;
* the label set is an ENUMERATED list, so it is checkable with no ontology file on the machine.

The panel is git-ignored (`data/cache`), so every test that needs it carries `needs_local_data` and
SKIPS BY NAME where the store is absent -- CI runs on a fresh checkout. The pure functions are
tested over STUB rows instead, which is behaviour and not existence: `enumerate_set`, `control`,
`design` and `clusters_needed` all take their input as an argument for exactly that reason.
"""

from __future__ import annotations

import pytest

from genomeos.attribution import astroargmax as aa

PANEL_HOW = "AlphaGenome's saved track metadata lives under the git-ignored data/cache store"


# ----------------------------------------------------------------------- the bands, and their order


def test_the_bands_are_ratios_so_no_control_value_could_have_fixed_them() -> None:
    assert aa.RATIO_ENRICHED_FLOOR == 2.0
    assert aa.RATIO_DEPLETED_CEILING == 0.8
    assert aa.EQUIVALENCE_BAND == 0.25
    assert aa.ABSOLUTE_LEVEL_FLOOR == 0.10
    # the equivalence band is symmetric in the ratio, which is what "relative" has to mean here
    assert round((1 / (1 + aa.EQUIVALENCE_BAND)) * (1 + aa.EQUIVALENCE_BAND), 12) == 1.0
    # and the reason an absolute band was rejected is argued, not asserted
    assert "91.7" in aa.WHY_RELATIVE_AND_NOT_ABSOLUTE
    assert "vacuous" in aa.WHY_RELATIVE_AND_NOT_ABSOLUTE


def test_the_enrichment_reading_needs_a_level_as_well_as_a_ratio() -> None:
    reading = aa.READINGS["enriched"]
    assert "AND" in reading
    assert str(aa.ABSOLUTE_LEVEL_FLOOR) in reading
    assert str(aa.RATIO_ENRICHED_FLOOR) in reading
    # the 2026-10-02 failure is named with its two numbers, so the rule carries its own reason
    assert "6.0%" in aa.WHY_A_RATIO_FLOOR_IS_NOT_ENOUGH
    assert "94%" in aa.WHY_A_RATIO_FLOOR_IS_NOT_ENOUGH
    assert "6.0%" in aa.EVERY_RATIO_CARRIES_ITS_LEVEL


def test_the_inconclusive_reading_says_the_data_cannot_tell_and_not_no_information() -> None:
    text = aa.READINGS["inconclusive"]
    assert "THE DATA CANNOT TELL" in text
    assert "Not 'no information'" in text


def test_the_instrument_assumptions_name_the_degenerate_bootstrap_and_its_route() -> None:
    d = aa.INSTRUMENT_ASSUMPTIONS["degenerate_bootstrap"]
    assert "DEGENERATE" in d and "decides nothing" in d
    assert "Clopper-Pearson" in d and "Wilson" in d and "n_eff" in d
    assert "identical_resample_share" in aa.INSTRUMENT_ASSUMPTIONS
    assert aa.INSTRUMENT_ASSUMPTIONS["assumed_independent"].startswith("chromosomes")
    # and what is NOT independent is listed rather than left to a reader
    assert any("TILING" in s for s in aa.INSTRUMENT_ASSUMPTIONS["known_not_independent"])


# --------------------------------------------------------- the label set is enumerated, not resolved


def test_the_registered_set_is_27_distinct_curies_over_four_prefixes() -> None:
    curies = [c for c, _ in aa.REGISTERED_SET]
    assert len(curies) == 27
    assert len(set(curies)) == 27
    assert sorted(set(c.split(":")[0] for c in curies)) == ["CL", "NTR", "UBERON"]


def test_every_excluded_borderline_term_is_disjoint_from_the_registered_set() -> None:
    members = set(c for c, _ in aa.REGISTERED_SET)
    for group in (aa.CNS_CELL_LINES, aa.NEURAL_CREST_TERMS):
        assert not (set(c for c, _ in group) & members)
    # the widest alternative is a strict superset and the narrowest a strict subset
    widest = set(c for c, _ in aa.ALTERNATIVE_SETS["CNS_AND_ALL_NEURAL_LINES"])
    narrow = set(c for c, _ in aa.ALTERNATIVE_SETS["BRAIN_ONLY"])
    assert members < widest
    assert narrow < members
    assert members - narrow == set(c for c, _ in aa.CNS_SPINAL_CORD)


def test_all_five_panel_prefixes_are_decided_by_the_enumeration_or_a_borderline_call() -> None:
    """CLO contributes nothing, and that is RECORDED rather than left as an omission."""
    decided = set(c for c, _ in aa.REGISTERED_SET)
    for group in (aa.CNS_CELL_LINES, aa.NEURAL_CREST_TERMS):
        decided |= set(c for c, _ in group)
    for call in aa.BORDERLINE_CALLS.values():
        decided |= set(t.split(" ")[0] for t in call["terms"] if ":" in t)
    assert sorted(set(c.split(":")[0] for c in decided)) == ["CL", "CLO", "EFO", "NTR", "UBERON"]
    assert "CLO:0007045" in decided
    assert aa.BORDERLINE_CALLS["clo_prefix"]["decision"] == "OUT"


def test_the_four_required_borderline_calls_are_each_decided_with_a_reason() -> None:
    for key in ("peripheral_nerve_tibial_and_sciatic", "pituitary_gland", "retina", "neuroblastoma_lines"):
        call = aa.BORDERLINE_CALLS[key]
        assert call["terms"] and call["decision"] and len(call["reason"]) > 120
    assert aa.BORDERLINE_CALLS["peripheral_nerve_tibial_and_sciatic"]["decision"] == "OUT"
    assert "DOES NOT EXIST" in aa.BORDERLINE_CALLS["retina"]["decision"]
    assert aa.BORDERLINE_CALLS["spinal_cord"]["decision"].startswith("IN")


def test_the_descendant_query_is_argued_against_with_the_measured_coverage_gap() -> None:
    text = aa.NOT_BY_DESCENDANT_RESOLUTION
    assert "cl-basic.obo" in text and "go-basic.obo" in text
    assert "0 tracked files" in text
    assert "146 of 667" in text and "21.89%" in text
    assert "UBERON:0000955" in text


# ---------------------------------------------------------------- the pure functions, over stub rows


def _rows(*spec: tuple[str, str, int]) -> list[dict[str, str]]:
    """Stub panel rows: (curie, biosample, how many), each row a distinct name."""
    out: list[dict[str, str]] = []
    for curie, bios, n in spec:
        for i in range(n):
            out.append(
                {
                    "name": f"{curie} src {bios} {i}",
                    "ontology_curie": curie,
                    "biosample_name": bios,
                    "data_source": "encode",
                    "output": "rna_seq",
                }
            )
    return out


def test_enumerate_set_counts_rows_per_term_and_names_a_term_the_panel_lacks() -> None:
    rows = _rows(("UBERON:0002037", "cerebellum", 5), ("EFO:0002067", "K562", 95))
    got = aa.enumerate_set((("UBERON:0002037", "cerebellum"), ("UBERON:9999999", "invented")), rows)
    assert got["terms"] == 2
    assert got["panel_rows"] == 5
    assert got["panel_row_share"] == 0.05
    assert got["absent_from_the_panel"] == ["UBERON:9999999"]
    assert "not a property of the genome" in got["panel_row_share_is_a_panel_share"]


def test_the_substring_crosscheck_reports_both_directions_of_its_own_error() -> None:
    rows = _rows(
        ("UBERON:0011907", "gastrocnemius medialis", 3),  # matches 'astro' inside 'gastro'
        ("UBERON:0002240", "spinal cord", 5),  # CNS, matches nothing
        ("UBERON:0002037", "cerebellum", 2),  # CNS and matches
    )
    got = aa.panel_substring_crosscheck(
        (("UBERON:0002240", "spinal cord"), ("UBERON:0002037", "cerebellum")), rows
    )
    assert got["matching_rows"] == 5
    assert got["label"].startswith("PANEL SHARE")
    assert got["disagrees_with_the_enumeration_on"] == 8
    assert [e["curie"] for e in got["matched_but_not_cns"]] == ["UBERON:0011907"]
    assert [e["curie"] for e in got["cns_but_not_matched"]] == ["UBERON:0002240"]


def test_the_control_joins_by_curie_and_ignores_a_label_that_merely_reads_like_one() -> None:
    cen = {
        "rules": 1000,
        "distinct_labels": 3,
        "per_cell": {
            "astrocyte": {"rules": 40, "ontology_term": "CL:0000127"},
            "brain": {"rules": 10, "ontology_term": "UBERON:0000955"},
            # a label whose NAME says brain but whose CURIE is not in the set: must not count
            "brain_microvascular_endothelial_cell": {"rules": 500, "ontology_term": "CL:9999999"},
            "unlabelled": {"rules": 450, "ontology_term": None},
        },
    }
    got = aa.control((("CL:0000127", "astrocyte"), ("UBERON:0000955", "brain")), cen)
    assert got["rules_with_a_label_in_the_set"] == 50
    assert got["rate"] == 0.05
    assert got["census_labels_matched_count"] == 2
    assert got["census_labels_carrying_no_ontology_term"] == 1
    assert got["curies_in_the_set_with_no_census_label"] == []
    # the forbidden control is named in the record the control itself carries
    assert "NEVER its share of TRACK NAMES" in got["never_a_track_share"]
    assert "2026-10-02" in got["never_a_track_share"]


def test_the_control_names_a_set_member_the_census_never_labelled() -> None:
    cen = {"rules": 100, "per_cell": {"a": {"rules": 1, "ontology_term": "CL:0000127"}}}
    got = aa.control((("CL:0000127", "astrocyte"), ("NTR:0000427", "neurosphere")), cen)
    assert got["curies_in_the_set_with_no_census_label"] == ["NTR:0000427"]
    assert got["rules_with_a_label_in_the_set"] == 1


def test_design_and_power_are_in_cluster_units_and_degrade_with_correlation() -> None:
    sizes = {"chr1": 100, "chr2": 100, "chr3": 100, "chr4": 100}
    d = aa.design(sizes)
    assert d["clusters"] == 4 and d["requests"] == 400
    assert d["effective_clusters_kish"] == 4.0  # equal sizes lose nothing
    uneven = aa.design({"chrA": 370, "chrB": 10, "chrC": 10, "chrD": 10})
    assert uneven["effective_clusters_kish"] < 2.0
    a = aa.clusters_needed(0.05, 0.04, 0.0, 50.0)
    b = aa.clusters_needed(0.05, 0.04, 0.05, 50.0)
    c = aa.clusters_needed(0.05, 0.04, 0.2, 50.0)
    assert a < b < c
    with pytest.raises(ValueError):
        aa.clusters_needed(0.04, 0.05, 0.0, 50.0)  # alternative below the null


def test_the_power_table_marks_an_underpowered_comparison_inconclusive_by_registration() -> None:
    out = aa.power_in_cluster_units(0.02, {"chr1": 137, "chr21": 5, "chrX": 7})
    a = out["comparison_a_vs_the_control"]
    assert a["null_the_lower_bound_must_clear"] == round(0.02 * aa.RATIO_ENRICHED_FLOOR, 8)
    assert a["absolute_level_floor_must_also_be_cleared"] == aa.ABSOLUTE_LEVEL_FLOOR
    assert "UNDERPOWERED BY REGISTRATION" in a["and_the_honest_part"]
    assert set(a["per_intracluster_correlation"]) == {"icc_0.0", "icc_0.05", "icc_0.2"}
    b = out["comparison_b_vs_the_K562_arm"]
    assert "NOT known to this registration" in b["clusters"]


# ------------------------------------------------------------------ the limitation that never softens


def test_training_exposure_is_registered_in_its_strongest_form() -> None:
    t = aa.TRAINING_EXPOSURE
    assert "RECALL OF SEEN DATA" in t
    assert "THIS DESIGN CANNOT DISTINGUISH THE TWO" in t
    assert "may not be shortened" in t
    assert "does NOT support" in t


def test_the_k562_arm_registers_the_axis_asymmetry_before_either_rate_exists() -> None:
    k = aa.K562_ARM
    assert "371" in k and "667" in k
    assert "not label-neutral" in k
    assert "collapsed pair" in k


def test_the_axis_is_not_a_population_and_the_code_constant_is_corrected_by_measurement() -> None:
    assert "NOT a denominator" in aa.AXIS_NOT_POPULATION
    assert "1,232 ANSWERS" in aa.AXIS_NOT_POPULATION
    m = aa.MEASURED_CORRECTION_TO_A_CODE_CONSTANT
    assert "316" in m and "371" in m and "conflated" in m


def test_the_file_says_plainly_that_it_authorises_nothing() -> None:
    assert "does not authorise the run" in aa.NOT_AN_AUTHORISATION
    assert "satisfies exactly one of them" in aa.NOT_AN_AUTHORISATION
    assert any("committed ALONE" in r for r in aa.REFUSALS)
    assert any("ABSENT from the payload" in r for r in aa.REFUSALS)


def test_the_disclosure_names_the_one_exposure_with_its_figure() -> None:
    d = aa.DISCLOSURE
    assert "neural_cell (CL:0002319) at 9,955 of 440,589 rules" in d
    assert "ALREADY FIXED" in d
    assert "no K562 argmax distribution of any kind" in d
    assert "Not one of the 1,232 requests has been sent." in d


# --------------------------------------------------- the panel itself: skipped by name where absent


@pytest.mark.needs_local_data(aa.PANEL_RELATIVE, how=PANEL_HOW)
def test_the_panel_is_the_file_declared_by_hash() -> None:
    from genomeos import manifest as mf

    digest, size, _ = mf.sha256_of(aa.PANEL)
    assert (digest, size) == (aa.PANEL_SHA256, aa.PANEL_BYTES)
    rows = aa.panel_rows()
    assert len(rows) == 667
    s = aa.panel_summary(rows)
    assert s["distinct_names"] == 371
    assert s["rows_that_collapse_under_name"] == 296
    assert s["rows_by_prefix"] == {"CL": 238, "CLO": 1, "EFO": 130, "NTR": 15, "UBERON": 283}
    assert s["distinct_ontology_terms"] == 285


@pytest.mark.needs_local_data(aa.PANEL_RELATIVE, how=PANEL_HOW)
def test_panel_rows_refuses_a_file_that_is_not_the_declared_one(tmp_path, monkeypatch) -> None:
    other = tmp_path / "not-the-panel.csv"
    other.write_text("name,ontology_curie,output\nx,UBERON:0000955,rna_seq\n", encoding="utf-8")
    monkeypatch.setattr(aa, "PANEL", other)
    with pytest.raises(ValueError, match="a different file is a different axis"):
        aa.panel_rows()


@pytest.mark.needs_local_data(aa.PANEL_RELATIVE, how=PANEL_HOW)
def test_the_enumerated_row_counts_are_what_the_registration_states() -> None:
    rows = aa.panel_rows()
    assert aa.enumerate_set(aa.REGISTERED_SET, rows)["panel_rows"] == 50
    counts = {
        name: aa.enumerate_set(terms, rows)["panel_rows"] for name, terms in aa.ALTERNATIVE_SETS.items()
    }
    assert counts == {
        "BRAIN_ONLY": 44,
        "CNS_AND_CNS_LINES": 58,
        "CNS_AND_ALL_NEURAL_LINES": 66,
    }
    # the substring cross-check is wrong in BOTH directions on 30 of the 667 rows
    cc = aa.panel_substring_crosscheck(aa.REGISTERED_SET, rows)
    assert cc["matching_rows"] == 40
    assert sum(e["panel_rows"] for e in cc["matched_but_not_cns"]) == 10
    assert sum(e["panel_rows"] for e in cc["cns_but_not_matched"]) == 20


# ----------------------------------------------- the committed registration: absences, and the order

THE_PAYLOAD_AND_NOT_THE_FILE = """
These test the PAYLOAD BUILDER and not the committed JSON, deliberately. The registration is
committed ALONE -- no code, no count, no test -- so a test that asked `must_be_committed` for it
would have to be red in the commit that adds this file and green only one commit later, and a lane
that leaves a red commit behind has broken the tree for every peer sharing it. The builder is what
produces the file, so the order of its keys and the ABSENCE of every run-fillable field are
properties of the builder; the file's identity is carried by the registration commit's own sha.
"""


def _payload() -> dict:
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "astroargmax_register", root / "scripts" / "astroargmax_register.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.payload()


#: Every key a RUN would fill. None of them may exist, not even as null: a present-and-empty field
#: reads as a slot waiting for a number, which is a preview of an answer that does not exist.
FIELDS_A_RUN_WOULD_FILL = (
    "observed_rate",
    "rate_on_the_1232",
    "answers_read",
    "argmax_labels",
    "interval",
    "ratio",
    "reading",
    "verdict",
    "k562_rate",
    "requests_sent",
    "elements_with_a_cns_argmax",
)


def test_the_registration_is_written_where_astroruns_clause_reads_it() -> None:
    """The file name is not a convention here: astrorun refuses the run until THIS path is in HEAD."""
    import importlib.util
    from pathlib import Path

    from genomeos.attribution import astrorun

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "astroargmax_register", root / "scripts" / "astroargmax_register.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert astrorun.ASTROARGMAX_REGISTRATION.as_posix() == f"data/results/{mod.RESULT}.json"


@pytest.mark.needs_local_data(aa.PANEL_RELATIVE, how=PANEL_HOW)
def test_the_payload_holds_no_field_a_run_would_fill() -> None:
    body = _payload()
    for field in FIELDS_A_RUN_WOULD_FILL:
        assert field not in body, f"{field} exists, so the registration can be read as a preview"
    assert body["alphagenome_requests"] == 0
    assert "no request is sent" in body["money"]
    assert "No astrocyte answer exists" in body["status"]


@pytest.mark.needs_local_data(aa.PANEL_RELATIVE, how=PANEL_HOW)
def test_the_payload_puts_the_bands_before_the_control() -> None:
    keys = list(_payload())
    for earlier, later in (
        ("bands", "label_set_rule"),
        ("label_set_rule", "label_set_CNS_PRIMARY"),
        ("label_set_CNS_PRIMARY", "control"),
        ("control", "disclosure"),
    ):
        assert keys.index(earlier) < keys.index(later), f"{earlier} must precede {later}"
    assert keys.index("panel_share_beside_the_control") > keys.index("control")
    assert keys.index("power_stated_in_advance") > keys.index("control")


@pytest.mark.needs_local_data(aa.PANEL_RELATIVE, how=PANEL_HOW)
def test_the_control_is_a_share_of_compiled_rules_and_not_of_tracks() -> None:
    body = _payload()
    ctl = body["control"]
    assert ctl["rules_total"] == 440589 == body["census"]["rules"]
    assert ctl["source"] == aa.CENSUS_RELATIVE
    assert 0 < ctl["rate"] < 1
    # the panel share sits beside it, LABELLED, and is a different number
    share = body["panel_share_beside_the_control"]
    assert share["label"].startswith("PANEL SHARE")
    assert share["panel_rows"] == 667
    assert share["panel_share"] != ctl["rate"]
    assert body["label_set_CNS_PRIMARY"]["panel_row_share"] != ctl["rate"]
    # and every alternative set gets its own control, computed the same way
    alt = body["control_for_each_alternative_set"]
    assert set(alt) == set(aa.ALTERNATIVE_SETS)
    assert all(a["rules_total"] == 440589 for a in alt.values())


@pytest.mark.needs_local_data(aa.PANEL_RELATIVE, how=PANEL_HOW)
def test_the_payload_declares_the_git_ignored_panel_by_hash() -> None:
    body = _payload()
    d = body["panel_declared_by_hash"]
    assert d["git_ignored"] is True
    assert d["sha256"] == aa.PANEL_SHA256
    assert d["bytes"] == aa.PANEL_BYTES
    assert "FRESH CHECKOUT CANNOT RECOMPUTE" in d["a_fresh_checkout_cannot_recompute_this"]
    panel = [e for e in body["result_manifest"]["inputs"] if e["path"] == aa.PANEL_RELATIVE]
    assert len(panel) == 1 and panel[0]["git_ignored"] is True


@pytest.mark.needs_local_data(aa.PANEL_RELATIVE, how=PANEL_HOW)
def test_the_payload_manifest_validates_and_claims_no_coordinate() -> None:
    from genomeos import manifest as mf

    m = _payload()["result_manifest"]
    assert mf.validate({**m, "code": {"git_sha": "x", "dirty": False}}) == []
    assert m["coordinates"].startswith("n/a:")
    assert m["assembly"].startswith("n/a:")
    assert m["parameters"]["ratio_enriched_floor"] == aa.RATIO_ENRICHED_FLOOR
    assert len(m["parameters"]["label_set_CNS_PRIMARY"]) == 27


# ====================================================== AMENDMENT 1: the unit is an ELEMENT, imported


def _row(gene: str, drop: float, rise: float, drop_t: str, rise_t: str) -> dict:
    return {
        "gene": gene,
        "max_drop_log2fc": drop,
        "max_rise_log2fc": rise,
        "max_drop_tissue": drop_t,
        "max_rise_tissue": rise_t,
    }


def test_element_label_is_one_label_per_element_from_the_winning_target_row() -> None:
    """Two gene rows, one element, ONE label: the target row's tissue and not the other's."""
    rows = [
        _row("WEAK", -0.2, 0.1, "K562", "HepG2"),
        _row("TARGET", -1.4, 0.3, "astrocyte", "brain"),
    ]
    assert aa.element_label(rows) == "astrocyte"
    # the non-target row's own argmax tissue is NOT the element's label
    assert aa.element_label([rows[0]]) == "K562"


def test_element_label_returns_none_where_no_gene_moves_by_min_effect() -> None:
    assert aa.min_effect() == 0.1
    assert aa.element_label([_row("FLAT", -0.01, 0.02, "K562", "HepG2")]) is None
    assert aa.element_label([]) is None


def test_min_effect_is_imported_and_no_literal_of_it_lives_in_this_lanes_module() -> None:
    """A copied literal and a re-implemented rule are the same defect."""
    import inspect

    from genomeos.predict import enhancer_target

    assert aa.min_effect() == enhancer_target.MIN_EFFECT
    src = inspect.getsource(aa)
    assert "MIN_EFFECT = " not in src, "the floor must be read from enhancer_target, never copied"


def test_element_label_DELEGATES_rather_than_agreeing_by_coincidence(monkeypatch) -> None:
    """Substitute predict_target and the label FOLLOWS it. A re-implementation would not move."""
    from genomeos.predict import enhancer_target

    monkeypatch.setattr(enhancer_target, "predict_target", lambda rows, **k: {"tissue": "planted"})
    assert aa.element_label([_row("ANY", -9.0, 0.0, "astrocyte", "brain")]) == "planted"
    monkeypatch.setattr(enhancer_target, "predict_target", lambda rows, **k: None)
    assert aa.element_label([_row("ANY", -9.0, 0.0, "astrocyte", "brain")]) is None


def test_element_label_uses_compiles_own_context_function(monkeypatch) -> None:
    """The label function is the compiler's, so a rule's label and an answer's cannot diverge."""
    from genomeos.attribution import compile as co

    monkeypatch.setattr(co, "context", lambda name: f"ctx<{name}>")
    assert aa.element_label([_row("T", -1.0, 0.0, "astrocyte", "brain")]) == "ctx<astrocyte>"


def test_the_module_holds_no_argmax_over_drop_and_rise_of_its_own() -> None:
    """The one place the unit is computed must not re-derive the target finder's comparison."""
    import inspect

    src = inspect.getsource(aa)
    for forbidden in ("max_drop_log2fc", "max_rise_log2fc", "max_drop_tissue", "max_rise_tissue"):
        assert forbidden not in src, f"{forbidden} appears, so the argmax is being re-implemented"


def test_label_population_reports_the_excluded_count_beside_the_denominator() -> None:
    got = aa.label_population(
        {
            "E1": [_row("A", -1.0, 0.0, "astrocyte", "brain")],
            "E2": [_row("B", -0.02, 0.01, "K562", "HepG2")],  # no gene moves by MIN_EFFECT
            "E3": [_row("C", 0.0, 2.0, "K562", "cerebellum")],
        }
    )
    assert got["elements_read"] == 3
    assert got["elements_with_a_target"] == 1 + 1
    assert got["elements_with_no_target"] == 1
    assert got["elements_with_no_target_named"] == ["E2"]
    assert got["elements_read"] == got["elements_with_a_target"] + got["elements_with_no_target"]
    assert got["denominator_is"] == "elements_with_a_target"
    assert got["labels_by_element"] == {"E1": "astrocyte", "E3": "cerebellum"}
    assert got["min_effect"] == 0.1
    assert "never silently dropped" in got["excluded_count_is_reported"]


def test_the_amendment_text_names_the_defect_the_parts_and_its_own_blindness() -> None:
    assert "THE UNIT UNDEFINED" in aa.AMENDMENT_1
    assert "compile.py` line 851" in aa.AMENDMENT_1
    assert "predict_target" in aa.AMENDMENT_1
    assert "not a wording fix" in aa.AMENDMENT_1.lower()
    assert "IMPORTED and never re-implemented" in aa.AMENDMENT_1_PART_1
    assert "D40" in aa.AMENDMENT_1_PART_1
    assert "REPORTED BESIDE THE RATE" in aa.AMENDMENT_1_PART_2
    assert "DECIDES NOTHING" in aa.AMENDMENT_1_PART_3
    assert "UNOBTAINABLE" in aa.AMENDMENT_1_IS_BLIND
    assert "2026-10-02" in aa.AMENDMENT_1_IS_BLIND
    assert "NOT rewritten" in aa.AMENDMENT_1_THE_REGISTRATION_IS_NOT_EDITED


def test_the_amendment_moves_no_band_and_names_the_blob_it_amends() -> None:
    assert aa.AMENDS == "astroargmax_registration"
    assert len(aa.AMENDS_BLOB_SHA256) == 64
    assert "byte-identical" in aa.AMENDMENT_1_WHAT_DOES_NOT_MOVE
    # the bands are the committed ones and this file does not redefine them
    assert (aa.RATIO_ENRICHED_FLOOR, aa.ABSOLUTE_LEVEL_FLOOR) == (2.0, 0.10)
    assert (aa.EQUIVALENCE_BAND, aa.RATIO_DEPLETED_CEILING) == (0.25, 0.8)


def _amendment_payload() -> dict:
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "astroargmax_amend_1", root / "scripts" / "astroargmax_amend_1.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.payload()


def test_the_amendment_payload_identifies_the_blob_it_amends_and_holds_the_three_parts() -> None:
    body = _amendment_payload()
    a = body["amended_registration"]
    assert a["commit"] == aa.AMENDS_COMMIT
    assert a["blob_sha256"] == aa.AMENDS_BLOB_SHA256
    assert a["bytes"] > 0
    for part in ("part_1_the_unit", "part_2_the_denominator", "part_3_the_secondary"):
        assert len(body[part]) > 200
    # the defect is EVIDENCED on the committed bytes, not asserted
    ev = body["the_evidence_counted_on_the_committed_bytes"]
    assert ev["gene row"] == 2
    assert ev["predict_target"] == 0 and ev["predicted_coding"] == 0 and ev["per element"] == 0
    # and the imported identities are recorded, so the fix rests on the code
    ids = body["the_label_path_in_the_code"]["imported_identities"]
    assert ids["predict_target"] == "genomeos.predict.enhancer_target.predict_target"
    assert ids["context"] == "genomeos.attribution.compile.context"
    assert body["the_label_path_in_the_code"]["min_effect"] == 0.1


def test_the_amendment_payload_is_blind_and_holds_no_field_a_run_would_fill() -> None:
    body = _amendment_payload()
    for field in FIELDS_A_RUN_WOULD_FILL:
        assert field not in body, f"{field} exists, so the amendment can be read as a preview"
    assert body["alphagenome_requests"] == 0
    assert body["network_requests"] == 0
    assert "0 of the 1,232 requests sent" in body["what_was_read_of_the_outcome_when_this_was_written"]
    assert "THE BLIND KIND" in body["which_kind_of_amendment_this_is"]
    assert body["bands_unchanged"]["ratio_enriched_floor"] == aa.RATIO_ENRICHED_FLOOR
    assert len(body["bands_unchanged"]["label_set_CNS_PRIMARY"]) == 27


def test_the_amendment_payload_manifest_validates_and_counts_the_exclusion() -> None:
    from genomeos import manifest as mf

    m = _amendment_payload()["result_manifest"]
    assert mf.validate({**m, "code": {"git_sha": "x", "dirty": False}}) == []
    assert m["parameters"]["denominator"] == "elements_with_a_target"
    assert m["parameters"]["excluded_and_counted"] == "elements_with_no_target"
    assert any("COUNTED and reported, never silently dropped" in e for e in m["exclusions"])
    assert set(m["partitions"]) == {"elements_with_a_target", "elements_with_no_target", "gene_rows"}
