# SPDX-License-Identifier: AGPL-3.0-or-later
"""The repression lane: the blinding of gate 1, the imported floors, and the registered refutation.

Every test here runs on synthetic records. No compiled program, no benchmark table and no measured
cache is read, so none of these tests can be made to pass by a count coming out one way.
"""

from __future__ import annotations

import pytest

from genomeos.attribution import cell2, fresh
from genomeos.attribution import context_evidence as ce
from genomeos.attribution import measured as ms
from genomeos.attribution import not_open_profile as nop
from genomeos.attribution import repress2 as rp


def rule(**kw):
    """One synthetic `not_open_profile.Rule`; only the fields a test varies are named."""
    base = dict(
        source=nop.SOURCE_PREDICTED,
        element="E1",
        chrom="chr21",
        start=1_000_000,
        end=1_000_200,
        cell="K562",
        gene="GENE1",
        element_class="enhancer",
        activity_axis=rp.REPRESSES,
        origin_axis="unique",
        effect=0.1,
        effect_band="weak",
        distance=5_000,
    )
    base.update(kw)
    return nop.Rule(**base)


def row(element="E1", pairs=()):
    return {"id": element, "measured": {"crispri": {"pairs": list(pairs)}}}


def pair(gene="GENE1", cell="K562", outcome=ms.DECREASE):
    return {"gene": gene, "cell": cell, "outcome": outcome, "regulated": outcome == ms.DECREASE}


class TestNothingIsSetHere:
    def test_both_floors_are_imported_not_chosen(self):
        assert rp.LOCUS_FLOOR == fresh.LOCUS_FLOOR == cell2.POOLED_LOCUS_FLOOR
        assert rp.POSITIVE_FLOOR == fresh.POSITIVE_FLOOR

    def test_the_locus_rule_and_span_are_cell2s_word_for_word(self):
        assert rp.LOCUS_RULE == cell2.INDEPENDENT_LOCUS_RULE
        assert rp.LOCUS_SPAN == cell2.INDEPENDENT_LOCUS_SPAN

    def test_the_bands_and_the_states_are_the_projects_own(self):
        assert rp.nop.EFFECT_BANDS == nop.EFFECT_BANDS
        assert rp.nop.STRONG_EFFECT == nop.STRONG_EFFECT
        for _name, _hb, _lb, hs, ls in rp.COMPARISONS:
            assert hs in (ce.STATE_OPEN, ce.STATE_NOT_OPEN)
            assert ls in (ce.STATE_OPEN, ce.STATE_NOT_OPEN)

    def test_the_retained_cells_are_the_sweeps_own_four(self):
        from genomeos.attribution.crispri import MODEL_CELLS

        assert rp.RETAINED_CELLS == MODEL_CELLS
        assert len(rp.RETAINED_CELLS) == 4


class TestGateOneIsBlind:
    def test_the_predicate_is_membership_in_the_two_significant_labels(self):
        assert rp.SIGNIFICANT_OUTCOMES == (ms.DECREASE, ms.INCREASE)
        assert rp.significant({"outcome": ms.DECREASE})
        assert rp.significant({"outcome": ms.INCREASE})
        for other in (ms.NULL_INFORMATIVE, ms.NULL_INCONCLUSIVE, ms.MISSING, None):
            assert not rp.significant({"outcome": other})

    def test_a_decrease_and_an_increase_are_counted_identically(self):
        """The whole point of the blinding: swapping the direction cannot move the count."""
        rules = [rule()]
        down = rp.eligible_links(rules, [row(pairs=[pair(outcome=ms.DECREASE)])], rp.REPRESSES, True)
        up = rp.eligible_links(rules, [row(pairs=[pair(outcome=ms.INCREASE)])], rp.REPRESSES, True)
        assert down == up
        assert len(down) == 1

    def test_no_effect_size_or_p_value_field_is_read(self):
        """A pair carrying a large effect size but a null outcome is not eligible."""
        p = pair(outcome=ms.NULL_INFORMATIVE)
        p["effect_size"] = -0.9
        p["p_adjusted"] = 0.0
        assert rp.eligible_links([rule()], [row(pairs=[p])], rp.REPRESSES, True) == []

    def test_the_gate_reads_only_the_blinded_fields(self):
        assert set(rp.GATE_FIELDS) == {
            "gene",
            "cell",
            "chrom",
            "start",
            "end",
            "outcome_is_one_of_the_two_significant_labels",
        }
        assert "effect_size" not in rp.GATE_FIELDS
        assert "p_adjusted" not in rp.GATE_FIELDS


class TestTheJoin:
    def test_a_different_gene_does_not_link(self):
        assert rp.eligible_links([rule()], [row(pairs=[pair(gene="OTHER")])], rp.REPRESSES, True) == []

    def test_a_different_element_does_not_link(self):
        assert rp.eligible_links([rule()], [row(element="E2", pairs=[pair()])], rp.REPRESSES, True) == []

    def test_a_different_cell_links_only_when_the_cell_is_ignored(self):
        rows = [row(pairs=[pair(cell="HepG2")])]
        assert rp.eligible_links([rule()], rows, rp.REPRESSES, True) == []
        assert len(rp.eligible_links([rule()], rows, rp.REPRESSES, False)) == 1

    def test_the_other_activity_axis_is_not_counted(self):
        rows = [row(pairs=[pair()])]
        assert rp.eligible_links([rule()], rows, rp.ACTIVATES, True) == []

    def test_a_combined_axis_is_in_neither_population(self):
        r = rule(activity_axis=f"{rp.ACTIVATES}, {rp.REPRESSES}")
        rows = [row(pairs=[pair()])]
        assert rp.eligible_links([r], rows, rp.REPRESSES, True) == []
        assert rp.eligible_links([r], rows, rp.ACTIVATES, True) == []

    def test_a_measured_rules_element_suffix_is_stripped_for_the_lookup(self):
        r = rule(element="E1_measured", source=nop.SOURCE_MEASURED, effect=None, effect_band=None)
        assert len(rp.eligible_links([r], [row(pairs=[pair()])], rp.REPRESSES, True)) == 1

    def test_the_key_carries_nothing_but_the_blinded_fields(self):
        keys = rp.eligible_links([rule()], [row(pairs=[pair()])], rp.REPRESSES, True)
        assert keys[0] == cell2.LocusKey("K562", "chr21", 1_000_000, 1_000_200, "GENE1")
        assert keys[0]._fields == cell2.LocusKey._fields


class TestTheGateReading:
    def test_below_either_floor_is_a_no_go(self):
        keys = [
            cell2.LocusKey("K562", "chr21", i * 5_000_000, i * 5_000_000 + 100, f"G{i}") for i in range(5)
        ]
        g = rp.gate(keys, "five links on five loci")
        assert g["independent_loci"] == 5
        assert not g["meets_both_floors"]
        assert g["reading"] == rp.GATE_NO_GO
        assert len(g["short_by"]) == 2

    def test_at_both_floors_exactly_is_a_pass(self):
        keys = [
            cell2.LocusKey("K562", "chr21", i * 5_000_000, i * 5_000_000 + 100, f"G{i}")
            for i in range(rp.LOCUS_FLOOR)
        ]
        keys += keys[: rp.POSITIVE_FLOOR - rp.LOCUS_FLOOR]
        g = rp.gate(keys, "at the floors")
        assert g["eligible_measured_links"] == rp.POSITIVE_FLOOR
        assert g["independent_loci"] == rp.LOCUS_FLOOR
        assert g["meets_both_floors"]
        assert g["reading"] == rp.GATE_PASS

    def test_links_within_the_span_are_one_locus(self):
        keys = [
            cell2.LocusKey("K562", "chr21", 1_000_000 + i * 1_000, 1_000_100 + i * 1_000, f"G{i}")
            for i in range(40)
        ]
        g = rp.gate(keys, "forty links inside one megabase")
        assert g["eligible_measured_links"] == 40
        assert g["independent_loci"] == 1
        assert g["reading"] == rp.GATE_NO_GO

    def test_the_two_readings_are_the_only_ones(self):
        assert len(rp.GATE_READINGS) == 2
        for word in ("validat", "confirm", "prove"):
            for r in rp.GATE_READINGS:
                assert word not in r.lower()


class TestThePredictionCanComeOutAgainstItself:
    @staticmethod
    def counts(shares):
        """(band, state, axis) counts giving each cell of the table the repression share asked for."""
        out = {}
        for key, share in shares.items():
            band, state = key.split("|")
            out[(band, state, rp.REPRESSES)] = int(round(share * 1000))
            out[(band, state, rp.ACTIVATES)] = 1000 - int(round(share * 1000))
        return out

    def rising(self):
        return self.counts(
            {
                f"strong|{ce.STATE_OPEN}": 0.10,
                f"strong|{ce.STATE_NOT_OPEN}": 0.20,
                f"weak|{ce.STATE_OPEN}": 0.30,
                f"weak|{ce.STATE_NOT_OPEN}": 0.40,
            }
        )

    def flat(self):
        return self.counts(
            {f"{b}|{s}": 0.25 for b in ("strong", "weak") for s in (ce.STATE_OPEN, ce.STATE_NOT_OPEN)}
        )

    def falling(self):
        return self.counts(
            {
                f"strong|{ce.STATE_OPEN}": 0.40,
                f"strong|{ce.STATE_NOT_OPEN}": 0.30,
                f"weak|{ce.STATE_OPEN}": 0.20,
                f"weak|{ce.STATE_NOT_OPEN}": 0.10,
            }
        )

    def read(self, counts):
        return rp.prediction_reading(rp.comparisons(rp.table(counts)))

    def test_a_rising_table_supports_it(self):
        r = self.read(self.rising())
        assert r["supported"]
        assert r["reading"] == rp.PREDICTION_SUPPORTED
        assert r["which_failed"] == []

    def test_a_flat_table_refutes_it(self):
        r = self.read(self.flat())
        assert not r["supported"]
        assert r["reading"] == rp.PREDICTION_REFUTED
        assert len(r["which_failed"]) == 4

    def test_a_falling_table_refutes_it(self):
        r = self.read(self.falling())
        assert not r["supported"]
        assert r["reading"] == rp.PREDICTION_REFUTED
        assert len(r["which_failed"]) == 4

    def test_one_failing_comparison_is_enough_to_refute_it(self):
        counts = self.counts(
            {
                f"strong|{ce.STATE_OPEN}": 0.10,
                f"strong|{ce.STATE_NOT_OPEN}": 0.05,  # falls where it should rise
                f"weak|{ce.STATE_OPEN}": 0.30,
                f"weak|{ce.STATE_NOT_OPEN}": 0.40,
            }
        )
        r = self.read(counts)
        assert not r["supported"]
        assert r["which_failed"] == ["not_open_above_open_within_strong"]

    def test_an_empty_cell_cannot_read_as_support(self):
        counts = self.counts({f"weak|{ce.STATE_NOT_OPEN}": 0.40})
        r = self.read(counts)
        assert not r["supported"]
        assert r["comparisons_on_an_empty_cell"]

    def test_every_share_carries_its_population(self):
        cells = rp.table(self.rising())
        assert len(cells) == 4
        for key, cell in cells.items():
            assert cell["assessable_rules"] == cell[rp.REPRESSES] + cell[rp.ACTIVATES]
            assert str(cell["assessable_rules"]) in cell["population_of_that_share"]
            assert cell["repression_share_ci95"] is not None, key

    def test_the_two_readings_are_the_only_ones_and_neither_says_validation(self):
        assert len(rp.PREDICTION_READINGS) == 2
        for r in rp.PREDICTION_READINGS:
            assert "validat" not in r.lower()


class TestTheInternalCheck:
    def test_a_call_agreeing_everywhere(self):
        signs = rp.cell_signs(0.5, {c: 0.1 for c in rp.RETAINED_CELLS})
        assert set(signs.values()) == {rp.AGREES}
        c = rp.concordance([signs])
        assert c["agrees_in_all_four"] == 1
        assert c["disagrees_in_all_four"] == 0

    def test_a_call_disagreeing_everywhere_is_counted_as_an_outlier(self):
        signs = rp.cell_signs(0.5, {c: -0.1 for c in rp.RETAINED_CELLS})
        assert set(signs.values()) == {rp.DISAGREES}
        c = rp.concordance([signs])
        assert c["disagrees_in_all_four"] == 1
        assert c["disagrees_in_all_four_share"] == 1.0
        assert "the 1 rules of this axis" in c["population_of_that_share"]

    def test_a_value_of_exactly_zero_has_no_sign(self):
        signs = rp.cell_signs(0.5, {"K562": 0.0, "HepG2": 0.1, "GM12878": -0.1, "IMR-90": 0.2})
        assert signs["K562"] == rp.NO_SIGN
        c = rp.concordance([signs])
        assert c["per_cell_value_exactly_zero"]["K562"] == 1
        assert c["disagrees_in_all_four"] == 0
        assert c["agrees_in_all_four"] == 0

    def test_an_element_missing_a_cell_is_not_scored_in_all_four(self):
        signs = rp.cell_signs(0.5, {"K562": 0.1})
        c = rp.concordance([signs])
        assert c["rules"] == 1
        assert c["rules_scored_in_all_four_retained_cells"] == 0
        assert c["rules_not_scored_in_all_four"] == 1
        assert c["disagrees_in_all_four_share"] is None

    def test_the_agreement_histogram_sums_to_the_scored_rules(self):
        rows = [
            rp.cell_signs(0.5, {"K562": 0.1, "HepG2": 0.1, "GM12878": -0.1, "IMR-90": -0.1}),
            rp.cell_signs(-0.5, {c: -0.1 for c in rp.RETAINED_CELLS}),
        ]
        c = rp.concordance(rows)
        assert sum(c["cells_agreeing_with_the_compiled_sign"].values()) == 2
        assert c["cells_agreeing_with_the_compiled_sign"]["2"] == 1
        assert c["cells_agreeing_with_the_compiled_sign"]["4"] == 1


@pytest.fixture
def reg():
    return rp.registration()


class TestTheRegistration:
    def test_it_says_it_is_not_a_result(self, reg):
        assert "registered before any count" in reg["status"]
        assert reg["alphagenome_requests"] == 0

    def test_the_mechanism_the_prediction_and_the_refutation_are_all_present(self, reg):
        assert reg["mechanism"] == rp.MECHANISM
        assert reg["prediction"] == rp.PREDICTION
        assert reg["refutation"] == rp.REFUTATION
        assert "refutes" in reg["refutation"]
        assert "does not rise" in reg["refutation"]

    def test_it_carries_the_limitations_word_for_word(self, reg):
        assert reg["limitations"]["not_closed"] == ce.NOT_CLOSED
        assert reg["limitations"]["not_validation"] == ce.NOT_VALIDATION

    def test_nothing_in_it_is_labelled_validation(self, reg):
        import json

        text = json.dumps(reg).lower()
        for sentence in text.split("."):
            if "validation" in sentence:
                assert "never" in sentence or "not " in sentence, sentence

    def test_the_measured_axis_caveat_is_on_the_record(self, reg):
        note = reg["limitations"]["measured_axis_is_not_the_link_direction"]
        assert "Regulated" in note
        assert "activates" in note

    def test_it_names_the_four_comparisons(self, reg):
        assert len(reg["comparisons"]) == 4
        assert len(set(reg["comparisons"])) == 4

    def test_the_two_suspended_results_carry_their_label(self, reg):
        assert set(reg["suspended_inputs"]) == {
            "data/results/context_evidence.json",
            "data/results/context_evidence_baserate.json",
        }
        for label in reg["suspended_inputs"].values():
            assert label == "suspended: not yet reproduced"
        assert "files_entry" in reg["suspension_call"]

    def test_it_says_no_figure_is_inherited_from_a_suspended_file(self, reg):
        note = reg["not_inherited_from_a_suspended_file"]
        assert "code and not a result" in note
        assert "state_for" in note

    def test_the_quoted_figures_name_their_file_and_scope(self, reg):
        q = reg["quoted_from"]
        assert q["file"] == "data/results/not_open_profile.json"
        assert "24 chromosomes" in q["scope"]
        assert "0.8511" in q["figures_quoted"]
