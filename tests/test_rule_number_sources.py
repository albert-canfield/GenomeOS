# SPDX-License-Identifier: AGPL-3.0-or-later
"""The rule-number source census of data/demo/gastrulation.bio.

These tests guard the shape of the census, not its answer: the population is read from the program,
the registration assigns nothing, and a finding may not be a bare verdict. The one test that does
constrain the answer is the one that refuses to let a fitted number be counted as measured.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.lang import rule_number_sources as rns

REGISTRATION = Path("data/results/rule_number_sources_registration.json")
CENSUS = Path("data/results/rule_number_sources_census.json")


def test_population_is_read_from_the_program_and_is_21():
    pop = rns.slots()
    assert len(pop) == 21, "3 fields x 7 rule declarations, enumerated from the file"
    assert len({(s["rule"], s["field"]) for s in pop}) == 21
    assert {s["field"] for s in pop} == set(rns.RULE_FIELDS)
    # every slot is written explicitly, so none of the 21 is a language default
    assert all(s["written_value"] is not None for s in pop)


def test_every_rule_of_the_program_writes_all_three_numbers():
    pop = rns.slots()
    by_rule: dict[str, set[str]] = {}
    for s in pop:
        by_rule.setdefault(s["rule"], set()).add(s["field"])
    assert len(by_rule) == 7
    assert all(fields == set(rns.RULE_FIELDS) for fields in by_rule.values())


def test_registration_assigns_no_class():
    payload = rns.registration()
    assert payload["findings_at_registration"] == 0
    assert payload["sources_fetched_at_registration"] == 0
    assert payload["verdict_at_registration"] == ""
    # the registration module's tables stay the tables AS REGISTERED: the filled ones live in
    # rule_number_sources_findings, and the writer loads them in its own process only
    assert rns.FINDINGS == {} and rns.SOURCES_FETCHED == () and rns.VERDICT == ""
    if REGISTRATION.exists():
        committed = json.loads(REGISTRATION.read_text())
        assert committed["findings_at_registration"] == 0
        assert committed["sources_fetched_at_registration"] == 0
        assert committed["verdict_at_registration"] == ""
        assert "verdict" not in committed and "counts" not in committed
    text = json.dumps(payload)
    for slot in payload["population"]:
        assert "class" not in slot
    # the class names appear only as definitions, never as an assignment on a slot
    assert all(isinstance(v, str) for v in payload["classes"].values())
    assert "measured_count" not in text


def test_the_cascade_is_the_classes_in_order():
    assert set(rns.CASCADE) == set(rns.CLASSES)
    assert len(rns.CASCADE) == len(rns.CLASSES)
    assert rns.CASCADE[0] == "M_measured_quoted" and rns.CASCADE[-1] == "U_unsourced"


def test_fitted_is_declared_never_measured_before_the_search():
    payload = rns.registration()
    assert "fitting" in payload["fitted_is_never_measured"]
    assert "never" in payload["fitted_is_never_measured"]
    # and the consequence for Hill coefficients is declared, not discovered
    assert "F_fitted_or_modelled" in rns.QUANTITY_DEFINITIONS["hill"]
    assert "FITTING" in rns.QUANTITY_DEFINITIONS["hill"]


def test_the_reported_header_line_is_still_false_and_still_unedited():
    """Line 3 of the program claims mutual repression; SOX17 is repressed by nothing."""
    text = Path(rns.PROGRAM).read_text()
    assert "Mutual repression makes the choice sharp" in text
    targets_of_sox17 = {s["rule"].split()[-1] for s in rns.slots() if s["rule"].startswith("Sox17 ")}
    assert targets_of_sox17 == {"TBXT", "SOX2"}
    repressors_of_sox17 = [s["rule"] for s in rns.slots() if s["rule"].endswith("inhibits SOX17")]
    assert repressors_of_sox17 == []
    assert any("Mutual repression" in r for r in rns.REPORTS)


def test_the_runtime_refuses_a_threshold_in_physical_units_for_this_program():
    """Axis 3 is not a stylistic point: a measured EC50 cannot be written into this program.

    The language parses a unit on a threshold; the unlocated GRN runtime refuses any rule that
    carries one. So for the seven thresholds, a measurement is not merely absent.
    """
    from genomeos.ir.model import Action
    from genomeos.lang.parser import parse_file
    from genomeos.runtime.grn import NetworkRuntime

    module = parse_file(rns.PROGRAM)
    rules = [r for r in module.rules if r.action in (Action.ACTIVATE, Action.INHIBIT)]
    assert rules and all(r.threshold_unit == "" for r in rules), "as written, every threshold is an amount"
    NetworkRuntime(module)  # runs as written
    rules[0].threshold_unit = "nM"
    with pytest.raises(ValueError, match="concentration"):
        NetworkRuntime(module)


@pytest.mark.skipif(not REGISTRATION.exists(), reason="registration not written yet")
def test_the_written_registration_matches_the_code_that_will_apply_it():
    d = json.loads(REGISTRATION.read_text())
    payload = rns.registration()
    for key in ("question", "classes", "cascade", "denominator", "population_size", "population"):
        assert d[key] == payload[key], key
    assert d["findings_at_registration"] == 0
    assert d["result_manifest"]["complete"] is True


@pytest.mark.skipif(not CENSUS.exists(), reason="census not written yet")
def test_every_slot_is_classified_and_every_finding_names_what_was_fetched():
    d = json.loads(CENSUS.read_text())
    findings = d["findings"]
    assert len(findings) == 21
    assert {(f["rule"], f["field"]) for f in findings} == {(s["rule"], s["field"]) for s in rns.slots()}
    for f in findings:
        for key in rns.REQUIRED_FINDING_KEYS:
            assert key in f, (f["rule"], f["field"], key)
        assert f["class"] in rns.CLASSES
        assert f["cited_source_attribution"] in rns.ATTRIBUTION
        assert f["commensurable"] in rns.COMMENSURABILITY
        assert f["source_does_not_claim"], "every entry says what its source does not claim"
        if f["class"] != "U_unsourced":
            assert f["fetched_url"], "a class above U rests on something fetched"
            assert f["fetched_what"]


@pytest.mark.skipif(not CENSUS.exists(), reason="census not written yet")
def test_the_counts_in_the_census_are_the_findings_recounted():
    d = json.loads(CENSUS.read_text())
    findings = d["findings"]
    counts: dict[str, int] = {c: 0 for c in rns.CLASSES}
    for f in findings:
        counts[f["class"]] += 1
    assert d["counts"] == counts
    assert sum(counts.values()) == 21 == d["population_size"]
    assert d["measured_strict"] == counts["M_measured_quoted"]
    assert d["measured_loose"] >= d["measured_strict"]


@pytest.mark.skipif(not CENSUS.exists(), reason="census not written yet")
def test_no_fitted_value_is_counted_as_measured():
    d = json.loads(CENSUS.read_text())
    fitted = [f for f in d["findings"] if f["class"] == "F_fitted_or_modelled"]
    assert len(fitted) == d["counts"]["F_fitted_or_modelled"]
    for f in fitted:
        assert f["fitted_how"], "a fitted entry says how the value was obtained"
    assert d["measured_strict"] + 0 == d["counts"]["M_measured_quoted"]
    # a Hill slot may never be measured: the class is declared by construction at registration
    assert not [f for f in d["findings"] if f["field"] == "hill" and f["class"] == "M_measured_quoted"]


@pytest.mark.skipif(not CENSUS.exists(), reason="census not written yet")
def test_the_census_changed_no_number_in_the_program():
    """The census records the written value of every slot; they are the program's own."""
    d = json.loads(CENSUS.read_text())
    written = {(f["rule"], f["field"]): f["written_value"] for f in d["findings"]}
    assert written == {(s["rule"], s["field"]): s["written_value"] for s in rns.slots()}
