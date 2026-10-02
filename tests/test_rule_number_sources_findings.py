# SPDX-License-Identifier: AGPL-3.0-or-later
"""The filled findings of the rule-number source census.

These tests guard what the findings may and may not do: every one rests on a URL that is in the
fetched list, the one attribution value added after registration is marked as added, the mappings
that were declined stay declined, and nothing fitted is counted as measured.
"""

from __future__ import annotations

import json
from pathlib import Path

from genomeos.lang import rule_number_sources as rns
from genomeos.lang import rule_number_sources_findings as rnf
from tests.committed_data import committed

REGISTRATION = Path("data/results/rule_number_sources_registration.json")
CENSUS_RELATIVE = "data/results/rule_number_sources_census.json"
CENSUS = Path(CENSUS_RELATIVE)


def test_the_filled_table_covers_every_slot_and_leaves_the_registered_one_empty():
    assert set(rnf.FINDINGS) == {(s["rule"], s["field"]) for s in rns.slots()}
    assert len(rnf.FINDINGS) == 21
    assert rns.FINDINGS == {}, "the registered table stays empty; this module holds the filled one"


def test_the_counts_are_recounted_from_the_findings():
    counts = rnf.counts()
    assert sum(counts.values()) == 21
    assert set(counts) == set(rns.CLASSES)
    assert counts["M_measured_quoted"] == 0, (
        "if this ever becomes non-zero the entry must quote a measured value from a fetched source, "
        "and a fitted or modelled value may never be what moved it"
    )


def test_no_hill_slot_is_measured_and_every_fitted_entry_says_how():
    for (rule, field), record in rnf.FINDINGS.items():
        if field == "hill":
            assert record["class"] != "M_measured_quoted", (rule, field)
        if record["class"] == "F_fitted_or_modelled":
            assert record["fitted_how"] and record["why_not_measured"], (rule, field)
            assert record["value_in_the_source"], (rule, field)


def test_every_finding_rests_on_a_url_in_the_fetched_list():
    fetched = {s["url"] for s in rnf.SOURCES_FETCHED}
    extra = {"https://pmc.ncbi.nlm.nih.gov/articles/PMC1885807/"}
    for (rule, field), record in rnf.FINDINGS.items():
        assert record["fetched_url"], (rule, field)
        assert set(record["fetched_url"]) <= fetched | extra, (rule, field)
        assert record["source_claims"] and record["source_does_not_claim"], (rule, field)
    for s in rnf.SOURCES_FETCHED:
        assert s["url"].startswith("https://"), s["key"]
        assert s["retrieved"] and s["claims"] and s["does_not_claim"], s["key"]


def test_the_added_attribution_value_is_marked_as_added_and_the_registration_lacks_it():
    added = rns.ATTRIBUTION_ADDED_AFTER_REGISTRATION
    assert set(added) == {"names_both_but_not_this_direction"}
    for name in added:
        assert rns.ATTRIBUTION[name].startswith("ADDED AFTER REGISTRATION")
    used = {record["cited_source_attribution"] for record in rnf.FINDINGS.values()}
    assert set(added) <= used, "a value is added only because a finding needed it"
    if REGISTRATION.exists():
        committed = json.loads(REGISTRATION.read_text())
        assert set(committed["cited_source_attribution_values"]) == set(rns.ATTRIBUTION) - set(added)
        assert len(committed["cited_source_attribution_values"]) == 5


def test_the_declined_mappings_are_recorded_and_used_for_nothing():
    """A measured number for a different pair is not a measured number for these."""
    assert len(rnf.DECLINED_MAPPINGS) == 2
    for m in rnf.DECLINED_MAPPINGS:
        assert m["why"] and m["values"] and m["declined_for"] and m["source"]
        assert any(s["key"] == m["source"] for s in rnf.SOURCES_FETCHED)
    thresholds = [r for (_, f), r in rnf.FINDINGS.items() if f == "threshold"]
    assert len(thresholds) == 7
    assert all(r["class"] in ("E_existence_only", "U_unsourced") for r in thresholds)
    assert all(r["commensurable"] == "arbitrary_units_no_conversion" for r in thresholds)


def test_the_one_rule_with_no_source_for_its_interaction_is_unsourced_on_all_three_numbers():
    unsourced = {rule for (rule, _), r in rnf.FINDINGS.items() if r["class"] == "U_unsourced"}
    assert unsourced == {"Sox17 inhibits SOX2"}
    for field in rns.RULE_FIELDS:
        record = rnf.FINDINGS[("Sox17 inhibits SOX2", field)]
        assert record["cited_source_attribution"] == "no_source_cited"


def test_the_verdict_states_the_count_and_claims_no_measurement():
    assert rnf.VERDICT.startswith("Of the 21 rule numbers")
    assert "0 have a measured value with a quoted source" in rnf.VERDICT
    assert "No number in the program was changed." in rnf.VERDICT


@committed(CENSUS_RELATIVE)
def test_the_written_census_matches_the_findings_module():
    d = json.loads(CENSUS.read_text())
    assert d["counts"] == rnf.counts()
    assert d["measured_strict"] == 0 and d["measured_loose"] == 0
    assert d["verdict"] == rnf.VERDICT
    assert len(d["sources_fetched"]) == len(rnf.SOURCES_FETCHED)
    assert len(d["attribution_findings"]) == 7
