# SPDX-License-Identifier: AGPL-3.0-or-later
"""The fresh non-K562 CRISPRi eligibility search: floors, blinding and the locus convention."""

from __future__ import annotations

import pytest

from genomeos.attribution import cell2, fresh


def candidate(**kwargs):
    base = dict(
        name="a set",
        cell_context="primary astrocytes",
        assay="CRISPRi",
        accession="doi x",
        pairs_tested=1000,
        positives=100,
        independent_loci=50,
        read_here="no outcome read here; established by grep over the record",
        access="open",
    )
    base.update(kwargs)
    return fresh.Candidate(**base)


class TestFloors:
    def test_the_locus_floor_is_cell2s_and_not_restated(self):
        assert fresh.LOCUS_FLOOR == cell2.POOLED_LOCUS_FLOOR == 20

    def test_the_positive_floor_is_thirty(self):
        assert fresh.POSITIVE_FLOOR == 30

    def test_the_excluded_cell_is_cell2s_primary(self):
        assert fresh.EXCLUDED_CELL == cell2.PRIMARY_CELL == "K562"

    def test_a_candidate_clearing_every_floor_is_eligible(self):
        assert candidate().eligible
        assert candidate().misses() == []

    def test_too_few_loci_misses_and_says_by_how_much(self):
        misses = candidate(independent_loci=12).misses()
        assert len(misses) == 1
        assert "independent loci 12" in misses[0] and "short of 20 by 8" in misses[0]

    def test_too_few_positives_misses_and_says_by_how_much(self):
        misses = candidate(positives=16).misses()
        assert len(misses) == 1
        assert "measured positives 16" in misses[0] and "short of 30 by 14" in misses[0]

    def test_a_locus_count_exactly_at_the_floor_clears_it(self):
        assert candidate(independent_loci=20, positives=30).eligible

    def test_k562_is_refused_however_large(self):
        misses = candidate(cell_context="K562", positives=10_000, independent_loci=500).misses()
        assert any("K562" in m for m in misses)

    def test_an_unknown_count_is_a_miss_not_a_pass(self):
        """A floor that cannot be checked has not been cleared."""
        assert not candidate(independent_loci=None).eligible
        assert "unchecked" in " ".join(candidate(independent_loci=None).misses())
        assert not candidate(positives=None).eligible
        assert "unchecked" in " ".join(candidate(positives=None).misses())

    def test_a_read_set_is_not_fresh_and_is_development_evidence(self):
        read = candidate(read_here="yes: its effect sizes were scored in 2026-09-27")
        assert not read.fresh
        assert not read.eligible
        assert any("not fresh" in m for m in read.misses())

    def test_freshness_needs_a_positive_statement(self):
        """Silence is not freshness: an unestablished record reads as not fresh."""
        assert not candidate(read_here="not established").fresh
        assert not candidate(read_here="").fresh
        assert not candidate(read_here="nothing found").fresh
        assert not candidate(read_here="no").fresh
        assert candidate(read_here="No outcome read here; established by grep").fresh


class TestBlinding:
    def test_permitted_columns_pass_through_unchanged(self):
        cols = ["Pair", "Hit", "EnhancerCoord", "GeneSymbol"]
        assert fresh.permit(cols) == cols

    @pytest.mark.parametrize(
        "column",
        [
            "log2FC",
            "Z_Sceptre",
            "P_Sceptre",
            "FDR_Sceptre",
            "GeneExpression",
            "Sensitivity_at_FC_0.15",
        ],
    )
    def test_an_outcome_magnitude_is_refused(self, column):
        with pytest.raises(fresh.BlindingError):
            fresh.permit(["Pair", column])

    def test_the_refusal_names_the_column_and_why(self):
        with pytest.raises(fresh.BlindingError) as exc:
            fresh.permit(["log2FC"])
        assert "log2FC" in str(exc.value)
        assert "outcome magnitude" in str(exc.value)

    def test_a_hit_call_is_permitted_because_counting_labels_is_the_point(self):
        assert fresh.permit(["Hit"]) == ["Hit"]

    def test_a_published_power_flag_is_permitted_but_a_sensitivity_figure_is_not(self):
        assert fresh.permit(["WellPowered_at_FC_0.15"]) == ["WellPowered_at_FC_0.15"]
        with pytest.raises(fresh.BlindingError):
            fresh.permit(["Sensitivity_at_FC_0.15"])

    def test_an_unknown_column_is_refused_even_without_a_forbidden_word(self):
        with pytest.raises(fresh.BlindingError):
            fresh.permit(["SomethingNew"])


class TestCoordinates:
    def test_a_well_formed_interval_parses(self):
        assert fresh.parse_coord("chr1:112390877-112391209") == ("chr1", 112390877, 112391209)

    def test_whitespace_is_tolerated(self):
        assert fresh.parse_coord("  chr21:100-200 ") == ("chr21", 100, 200)

    @pytest.mark.parametrize("text", ["chr1:100", "1:100-200", "chr1-100-200", "", "chr1:abc-200"])
    def test_anything_unspellable_raises_rather_than_being_skipped(self, text):
        with pytest.raises(ValueError):
            fresh.parse_coord(text)

    def test_a_reversed_interval_raises(self):
        with pytest.raises(ValueError):
            fresh.parse_coord("chr1:200-100")


class TestLocusGrouping:
    def test_the_grouping_is_cell2s(self):
        rows = [
            {"c": "chr1:1000-1100", "g": "AAA"},
            {"c": "chr1:1200-1300", "g": "BBB"},
            {"c": "chr5:9000-9100", "g": "CCC"},
        ]
        keys = fresh.locus_keys("astrocyte", rows, "c", "g")
        # the two chr1 elements are 100 bp apart, so they chain into one locus
        assert cell2.count_loci(keys) == 2

    def test_a_shared_measured_gene_is_one_locus_across_chromosomes(self):
        rows = [
            {"c": "chr1:1000-1100", "g": "SAME"},
            {"c": "chr9:5000000-5000100", "g": "SAME"},
        ]
        keys = fresh.locus_keys("astrocyte", rows, "c", "g")
        assert cell2.count_loci(keys) == 1

    def test_elements_beyond_the_span_are_separate_loci(self):
        far = cell2.INDEPENDENT_LOCUS_SPAN + 10_000
        rows = [
            {"c": "chr1:1000-1100", "g": "AAA"},
            {"c": f"chr1:{1000 + far}-{1100 + far}", "g": "BBB"},
        ]
        keys = fresh.locus_keys("astrocyte", rows, "c", "g")
        assert cell2.count_loci(keys) == 2

    def test_every_row_is_keyed_and_none_is_dropped(self):
        rows = [{"c": f"chr2:{i * 10_000_000}-{i * 10_000_000 + 100}", "g": f"G{i}"} for i in range(7)]
        assert len(fresh.locus_keys("astrocyte", rows, "c", "g")) == 7

    def test_one_unparseable_row_fails_the_whole_count(self):
        rows = [{"c": "chr1:1000-1100", "g": "AAA"}, {"c": "nonsense", "g": "BBB"}]
        with pytest.raises(ValueError):
            fresh.locus_keys("astrocyte", rows, "c", "g")


class TestVerdict:
    def test_a_verdict_with_one_eligible_set_is_a_go(self):
        verdict = fresh.Verdict([candidate(), candidate(name="small", positives=4)])
        out = verdict.to_dict()
        assert out["verdict"] == "go"
        assert out["eligible_sets"] == ["a set"]
        assert out["candidates_examined"] == 2
        assert out["eligible_count"] == 1

    def test_a_verdict_with_nothing_eligible_is_a_no_go(self):
        verdict = fresh.Verdict([candidate(positives=4), candidate(name="k", cell_context="K562")])
        out = verdict.to_dict()
        assert out["verdict"] == "no-go"
        assert out["eligible_sets"] == []

    def test_an_empty_search_is_a_no_go_rather_than_a_go(self):
        assert fresh.Verdict([]).to_dict()["verdict"] == "no-go"

    def test_the_verdict_carries_the_floors_and_the_registered_rule_verbatim(self):
        floors = fresh.Verdict([]).to_dict()["floors"]
        assert floors["independent_loci_at_least"] == 20
        assert floors["measured_positives_at_least"] == 30
        assert floors["cell_context_not"] == "K562"
        assert floors["locus_rule"] == cell2.INDEPENDENT_LOCUS_RULE
        assert floors["locus_span_bp"] == 1_000_000

    def test_the_floors_say_they_were_fixed_before_any_candidate_was_opened(self):
        assert fresh.Verdict([]).to_dict()["floors"]["fixed_before_any_candidate_was_opened"]

    def test_the_rule_is_not_sold_as_biological_independence(self):
        note = fresh.Verdict([]).to_dict()["floors"]["note"]
        assert "not established biological independence" in note

    def test_every_candidate_carries_how_freshness_was_established(self):
        row = candidate().to_dict()
        assert "grep" in row["read_here"]
        assert row["fresh"] is True
        assert row["misses"] == []
