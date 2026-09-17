# SPDX-License-Identifier: AGPL-3.0-or-later
"""The experimental layer: the overlap rule, the measured negatives and the census that counts them.

Everything here is synthetic and local. Nothing is fetched, and the assay caches are never read: the
`Layer` is built from hand-made rows so the rule can be checked at its edges.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import measured
from genomeos.attribution.compile import compile_chromosome
from genomeos.attribution.mpra import Element as MpraElement
from genomeos.attribution.vista import VistaElement
from genomeos.bio import evaluate, parse_tests
from genomeos.lang import parse

CHROM = "chrT"


def pair(start, end, gene, regulated, effect=-0.3, cell="K562"):
    return measured.CrispriPair(
        chrom=CHROM,
        start=start,
        end=end,
        gene=gene,
        cell=cell,
        dataset="Synthetic2026",
        reference="synthetic",
        regulated=regulated,
        significant=regulated,
        effect_size=effect,
        p_adjusted=0.001 if regulated else 0.9,
    )


def mpra_element(start, end, name, **activity):
    e = MpraElement(CHROM, start, end, name)
    e.activity.update(activity)
    return e


def satmut_element(start, end, experiment, functional=(), strong=(), repeats=1, bases_from=None):
    """A synthetic saturation-mutagenesis experiment: one base table over (start, end), 0-based.

    `bases_from` measures only part of the span, which is how the "matched but none of this element's
    own bases were measured" case is built - a real experiment's flanks can fall outside the element.
    """
    functional, strong = set(functional), set(strong)
    table = {
        p + 1: {
            "ref": "A",
            "measured": 3,
            "functional": p in functional,
            "strong": p in strong,
            "effect": 0.9 if p in strong else (0.1 if p in functional else 0.0),
        }
        for p in range(bases_from if bases_from is not None else start, end)
    }
    return measured.SatmutElement(CHROM, start, end, experiment, repeats, table)


# --- the rule ---------------------------------------------------------------------------------
def test_reciprocal_overlap_is_the_smaller_of_the_two_fractions():
    assert measured.reciprocal_overlap(1000, 1300, 1000, 1300) == 1.0
    assert measured.reciprocal_overlap(1000, 1300, 1300, 1600) == 0.0
    assert measured.reciprocal_overlap(1000, 1300, 900, 1400) == pytest.approx(300 / 500)
    # a long interval that swallows a short element scores low, because it is not a measurement of it
    assert measured.reciprocal_overlap(1000, 1300, 0, 10_000) == pytest.approx(300 / 10_000)


def test_containment_alone_never_meets_the_rule():
    # a 2 kb VISTA sequence around a 300 bp element: contained, but not measured
    assert measured.contains(1000, 1300, 500, 2500)
    assert not measured.measures(1000, 1300, 500, 2500, measured.RECIPROCAL_OVERLAP)
    assert not measured.measures(1000, 1300, 500, 2500, 0.25)


def test_the_rule_constant_is_the_one_the_module_uses():
    assert measured.RECIPROCAL_OVERLAP in measured.OVERLAP_SENSITIVITY
    assert measured.measures(1000, 1300, 1050, 1350)  # 250/300 both ways
    assert not measured.measures(1000, 1300, 1250, 1550)  # 50/300 both ways


# --- the layer --------------------------------------------------------------------------------
@pytest.fixture
def layer():
    return measured.Layer(
        chrom=CHROM,
        crispri=[
            pair(1000, 1300, "AAA", True, effect=-0.42),
            pair(1000, 1300, "BBB", False, effect=-0.01),
            pair(5000, 5300, "CCC", False, effect=0.02),
        ],
        lentimpra=[
            mpra_element(1010, 1290, "syn_1", K562=2.4, HepG2=0.1),
            mpra_element(5020, 5280, "syn_2", K562=-0.8, HepG2=-0.3),
        ],
        vista=[VistaElement("hs9001", CHROM, 9000, 9350, "positive", ("fb", "ht"), 1)],
        crispri_invalid=7,
    )


@pytest.fixture
def elements():
    """Four compiled elements: measured and regulated, measured and negative, contained only, blind."""

    def el(eid, start, end, gene, log2, action="activates"):
        return {
            "id": eid,
            "start": start,
            "end": end,
            "domain": "D1",
            "predicted_coding": {"gene": gene, "log2_fold_change": log2, "action": action},
        }

    return [
        el("E1", 1000, 1300, "AAA", -0.5),
        el("E2", 5000, 5300, "CCC", -0.4),
        el("E3", 9100, 9200, "DDD", 0.3, "represses"),
        el("E4", 80_000, 80_300, "EEE", -0.2),
    ]


def test_for_element_returns_only_measurements_of_that_element(layer, elements):
    got = layer.for_element(1000, 1300)
    assert set(got) == {"crispri", "lentimpra"}
    assert got["crispri"]["genes_regulated"] == ["AAA"]
    assert got["crispri"]["genes_not_regulated"] == ["BBB"]
    assert got["lentimpra"]["cells_active"] == ["K562"]
    assert layer.for_element(80_000, 80_300) == {}
    # E3 sits inside the VISTA element but is not measured by it
    assert layer.for_element(9100, 9200) == {}
    assert layer.for_element(9000, 9350)["vista"]["positive"] == ["hs9001"]


def test_a_measured_negative_is_a_row_and_keeps_its_evidence_kind(layer, elements):
    rows = measured.rows(CHROM, elements, layer)
    by_id = {r["id"]: r for r in rows}
    assert set(by_id) == {"E1", "E2"}  # E3 is only contained, E4 is untouched
    e2 = by_id["E2"]
    assert e2["measured"]["crispri"]["genes_regulated"] == []
    assert e2["measured"]["crispri"]["genes_not_regulated"] == ["CCC"]
    # the screen tested exactly the predicted gene and found nothing: that is a disagreement, not a gap
    assert e2["agreement"]["crispri"] == measured.DISAGREES
    assert e2["confidence"] == measured.PERTURBATION_CONFIDENCE
    # and it makes no rule, because a rule would state a relation the assay says is absent
    assert measured.rule_lines(e2) == []


def test_agreement_is_three_valued_and_names_the_untested_case(layer):
    m = layer.for_element(1000, 1300)
    assert measured.agreement("AAA", m)["crispri"] == measured.AGREES
    assert measured.agreement("BBB", m)["crispri"] == measured.DISAGREES
    assert measured.agreement("ZZZ", m)["crispri"] == measured.NOT_TESTED
    assert measured.agreement("ZZZ", m)["crispri_detail"] == "another gene regulated"


def test_the_two_agreement_rates_have_different_denominators(layer, elements):
    rows = measured.rows(CHROM, elements, layer)
    c = measured.census(CHROM, rows, elements, layer)
    assert c["agreement_denominator_all_matched_elements"] == 2
    assert c["agreement_denominator_predicted_gene_tested"] == 2
    assert c["agrees"]["crispri"] == 1
    assert c["disagrees"]["crispri"] == 1
    # the narrower rate is reported under its own name, never in place of the wider one
    assert c["agreement_rate_over_all_matched_elements"] == 0.5
    assert c["agreement_rate_where_the_predicted_gene_was_tested"] == 0.5


def test_eligibility_is_reported_before_the_census_and_is_a_weaker_bar(layer, elements):
    """What could have been found at all, so a small census reads as coverage and not as biology."""
    e = measured.eligibility(elements, layer)
    assert e["compiled_elements"] == 4
    # E1 and E2 are measured, E3 sits inside a VISTA test it does not meet the rule with: three of four
    assert e["elements_in_an_assay_footprint"] == 3
    assert e["elements_no_assay_ever_covered"] == 1  # E4, and it is in its own bucket
    assert e["share_eligible"] == 0.75
    assert e["by_assay"]["vista"] == 1
    assert e["rule"] == measured.ELIGIBILITY_RULE
    # eligibility can never be narrower than the overlap rule
    rows = measured.rows(CHROM, elements, layer)
    assert e["elements_in_an_assay_footprint"] >= len(rows)


def test_the_census_counts_coverage_negatives_and_two_kinds_of_fact(layer, elements):
    rows = measured.rows(CHROM, elements, layer)
    c = measured.census(CHROM, rows, elements, layer)
    assert c["compiled_elements"] == 4
    assert c["elements_in_an_assay_footprint"] == 3
    assert c["elements_no_assay_ever_covered"] == 1
    assert c["elements_with_any_measurement"] == 2
    # raised out of the eligible, not out of everything: two different denominators, two names
    assert c["raised_share_of_the_eligible"] == round(2 / 3, 5)
    assert c["coverage_elements_measured"] == 0.5
    # two names, because they answer two questions: what sits beside a prediction, and what was added
    assert c["facts_with_a_measured_counterpart"] == 4  # an element block and its rule, per element
    assert c["experimental_facts_added"] == 3  # two measured element blocks and one measured rule
    assert c["experimental_element_blocks"] == 2
    assert c["experimental_rule_blocks"] == 1
    assert c["crispri_pairs_measured_as_not_regulated"] == 2
    assert c["crispri_pairs_the_benchmark_calls_invalid"] == 7


def test_the_pooled_census_sums_what_can_be_summed(layer, elements):
    rows = measured.rows(CHROM, elements, layer)
    c = measured.census(CHROM, rows, elements, layer)
    pooled = measured.pool([c, c])
    assert pooled["chromosomes"] == 2
    assert pooled["compiled_elements"] == 8
    assert pooled["elements_with_any_measurement"] == 4
    assert pooled["coverage_elements_measured"] == 0.5
    assert pooled["experimental_rule_blocks"] == 2
    assert pooled["eligibility"]["elements_in_an_assay_footprint"] == 6
    assert pooled["eligibility"]["elements_no_assay_ever_covered"] == 2
    assert pooled["raised_share_of_the_eligible"] == round(4 / 6, 5)


def test_sensitivity_reports_the_rule_at_three_values_and_the_refused_containments(layer, elements):
    s = measured.sensitivity(elements, layer)
    assert s["used"] == measured.RECIPROCAL_OVERLAP
    assert list(s["at"]) == ["0.25", "0.5", "0.75"]
    counts = [s["at"][k]["any"] for k in ("0.25", "0.5", "0.75")]
    assert counts == sorted(counts, reverse=True)  # a stricter rule never matches more
    assert s["containment_not_upgraded"]["vista"] == 1  # E3 inside hs9001, and left predicted


# --- the base-level assay: which bases matter, which is not which elements do something ----------
@pytest.fixture
def base_elements():
    """Four elements against saturation mutagenesis: functional, inert, in the footprint, unmeasured."""

    def el(eid, start, end, gene):
        return {
            "id": eid,
            "start": start,
            "end": end,
            "domain": "D1",
            "predicted_coding": {"gene": gene, "log2_fold_change": -0.5, "action": "activates"},
        }

    return [
        el("S1", 1000, 1300, "AAA"),  # measured, and 50 of its bases matter
        el("S2", 5000, 5300, "CCC"),  # measured base by base, every one of them inert
        el("S3", 9100, 9200, "DDD"),  # inside the footprint, below the overlap rule
        el("S4", 6000, 6150, "EEE"),  # meets the rule, but its own bases were never measured
        el("S5", 80_000, 80_300, "FFF"),  # no assay ever looked here
    ]


@pytest.fixture
def base_layer():
    return measured.Layer(
        chrom=CHROM,
        satmut=[
            satmut_element(1000, 1300, "SYN1", functional=range(1100, 1150), strong=range(1100, 1110)),
            satmut_element(5000, 5300, "SYN2"),
            satmut_element(9180, 9480, "SYN3", functional=range(9180, 9280)),
            satmut_element(6000, 6300, "SYN4", functional=range(6200, 6300), bases_from=6200),
        ],
    )


def test_a_base_level_assay_can_agree_and_cannot_disagree(base_layer, base_elements):
    rows = measured.rows(CHROM, base_elements, base_layer)
    by_id = {r["id"]: r for r in rows}
    assert set(by_id) == {"S1", "S2", "S4"}  # S3 is only in the footprint, S5 was never looked at
    assert by_id["S1"]["agreement"]["satmut"] == measured.AGREES
    # every measured base inert: a measurement of those bases, and not a verdict on the element
    assert by_id["S2"]["agreement"]["satmut"] == measured.BASES_INERT
    assert by_id["S2"]["agreement"]["assays_disagreeing"] == 0
    assert by_id["S2"]["agreement"]["verdict"] == measured.NOT_TESTED
    # matched, but the experiment's measured bases lie outside this element: the never-looked zero
    assert by_id["S4"]["agreement"]["satmut"] == measured.BASES_NOT_MEASURED
    assert by_id["S4"]["measured"]["satmut"]["bases_measured"] == 0
    assert "measured none of its bases" in measured.basis_text(by_id["S4"])
    # and it is a reporter, so it carries the reporter's confidence and no rule
    assert by_id["S1"]["confidence"] == measured.REPORTER_CONFIDENCE
    assert measured.rule_lines(by_id["S1"]) == []


def test_what_is_said_is_said_over_the_elements_own_bases(base_layer, base_elements):
    s = base_layer.for_element(1000, 1300)["satmut"]
    assert (s["bases_in_element"], s["bases_measured"]) == (300, 300)
    assert (s["bases_functional"], s["bases_strong"]) == (50, 10)
    assert s["share_functional_of_measured"] == round(50 / 300, 4)
    assert s["share_of_the_element_measured"] == 1.0
    # a shorter element inside the same experiment is told only about its own bases
    half = base_layer.for_element(1000, 1150)["satmut"]
    assert (half["bases_measured"], half["bases_functional"]) == (150, 50)


def test_the_census_names_the_inert_and_the_never_raised_apart(base_layer, base_elements):
    rows = measured.rows(CHROM, base_elements, base_layer)
    c = measured.census(CHROM, rows, base_elements, base_layer)
    s = c["satmut"]
    assert c["eligibility"]["by_assay"]["satmut"] == 4  # S1, S2, S3, S4 - S5 is in nobody's footprint
    assert s["elements_measured_base_by_base"] == 3
    assert s["elements_with_a_functional_base"] == 1
    assert s["elements_measured_and_every_base_inert"] == 1  # measured and absent
    assert s["elements_matched_but_no_base_of_theirs_measured"] == 1  # never looked, same element
    assert s["elements_in_the_footprint_never_raised"] == 1  # S3, never looked, whole element
    assert (s["bases_measured"], s["bases_functional"], s["bases_strong"]) == (600, 50, 10)
    # the assay's zero in `disagrees` is a property of the assay, and the census says so in words
    assert c["disagrees"]["satmut"] == 0
    assert s["why_it_can_never_disagree"] == measured.SATMUT_CANNOT_DISAGREE
    assert measured.pool([c, c])["satmut"]["bases_functional"] == 100


def test_a_missing_satmut_cache_is_no_measurement_and_never_a_fetch(tmp_path):
    """The layer must be readable where no satmut cache is, and must then measure nothing."""
    assert measured.load_satmut(CHROM, tmp_path / "absent" / "elements.tsv.gz") == []


def test_the_program_states_the_base_level_fact_and_raises_no_rule(tmp_path, base_layer, base_elements):
    text = compile_chromosome(CHROM, _results_dir(tmp_path, base_elements), layer=base_layer)
    module = parse(text)
    assert module.entities["S1_measured"].evidence.kind == "experimental"
    # the inert element keeps `experimental` too: no effect on these bases is a measurement
    assert module.entities["S2_measured"].evidence.kind == "experimental"
    assert "no base measured here matters" in text
    assert not [r for r in module.rules if r.evidence.kind == "experimental"]
    # the header states the base-level count, and says which bucket the inert element is not in
    assert "neither the agreements nor the disagreements" in text
    assert all(r["ok"] for r in evaluate(module, parse_tests(text)))


# --- the compiled program ----------------------------------------------------------------------
def _results_dir(tmp_path: Path, elements: list[dict]) -> Path:
    d = tmp_path / "results"
    d.mkdir()
    (d / f"budget_{CHROM}.json").write_text(
        json.dumps(
            {
                "unknown_bp": 100_000,
                "constrained_fraction": 0.02,
                "blocks": [
                    {
                        "start": 0,
                        "end": 100_000,
                        "class": "regulatory",
                        "guess": {"tier": "constrained_unknown", "label": "real unknown", "confidence": 0.3},
                        "phylop": {"fraction_above": 0.4, "bases": 100_000},
                        "elements": {"n": 3},
                    }
                ],
            }
        )
    )
    (d / f"enhancer_targets_{CHROM}.json").write_text(json.dumps({"elements": elements}))
    return d


def test_the_program_carries_the_measurement_beside_the_prediction(tmp_path, layer, elements):
    text = compile_chromosome(CHROM, _results_dir(tmp_path, elements), layer=layer)
    module = parse(text)
    # two names, not one: the prediction and the measurement are separate blocks
    assert "E1" in module.entities and "E1_measured" in module.entities
    assert module.entities["E1"].evidence.kind == "predicted"
    assert module.entities["E1_measured"].evidence.kind == "experimental"
    assert module.entities["E1"].locus == module.entities["E1_measured"].locus
    # the measured negative is stated, with its evidence kind, and names the gene it found nothing on
    assert module.entities["E2_measured"].evidence.kind == "experimental"
    assert "measured no effect on CCC" in text
    # one experimental rule, for the one link a screen measured as regulated
    experimental = [r for r in module.rules if r.evidence.kind == "experimental"]
    assert len(experimental) == 1
    assert (experimental[0].source, experimental[0].action, experimental[0].target) == (
        "E1_measured",
        "activates",
        "AAA",
    )


def test_the_programs_own_test_lines_still_pass_with_the_layer(tmp_path, layer, elements):
    text = compile_chromosome(CHROM, _results_dir(tmp_path, elements), layer=layer)
    results = evaluate(parse(text), parse_tests(text))
    assert results and all(r["ok"] for r in results), results


def test_no_measurement_states_an_empty_layer_rather_than_staying_silent(tmp_path, elements):
    empty = measured.Layer(chrom=CHROM)
    assert empty.empty
    text = compile_chromosome(CHROM, _results_dir(tmp_path, elements), layer=empty)
    module = parse(text)
    assert not [e for e in module.entities if e.endswith("_measured")]
    assert not [r for r in module.rules if r.evidence.kind == "experimental"]
    # and it says so, so "nothing measured" cannot be confused with "the layer was never written"
    assert f"the experimental layer (0 of {len(elements)} elements" in text
    assert "states no experimental fact" in text
    assert all(r["ok"] for r in evaluate(module, parse_tests(text)))


# --- the census must not be able to report a zero it never read ---------------------------------
def test_a_census_that_reads_nothing_is_not_a_census_that_measured_nothing(tmp_path):
    """`evidence.collect` skips the compiled tree before parsing it, so the flag is load-bearing."""
    from genomeos.evidence import COMPILED_DIR, collect

    root = tmp_path / "root"
    (root / COMPILED_DIR).mkdir(parents=True)
    (root / COMPILED_DIR / f"noncoding_{CHROM}.bio").write_text(
        f"module human.noncoding.{CHROM}\n"
        'element X_measured { class: enhancer; locus: chrT:1-2; evidence: experimental "a screen"; '
        "confidence: 0.90 }\n"
    )
    blind = collect(root)
    assert blind["whole"]["facts"] == 0  # read nothing, and says so as an empty census
    seeing = collect(root, compiled=True)
    assert seeing["whole"]["facts"] == 1
    assert seeing["whole"]["by_evidence"]["experimental"] == 1


def test_the_committed_chromosome_program_states_experimental_facts():
    """chr21 is compiled from caches that are on this machine, so its layer must not be empty.

    If this fails after a recompile, the layer was dropped rather than measured away: check that the
    assay caches under data/knowledge are present before believing a zero.
    """
    path = Path("data/organisms/human/noncoding_chr21.bio")
    if not path.exists():
        pytest.skip("chr21 program not compiled on this machine")
    module = parse(path.read_text())
    experimental = [
        e
        for e in module.entities.values()
        if getattr(e, "evidence", None) and e.evidence.kind == "experimental"
    ]
    assert len(experimental) > 0, "the compiled chr21 program states no experimental fact"
    assert all(e.id.endswith("_measured") for e in experimental)
    assert all(r["ok"] for r in evaluate(module, parse_tests(path.read_text())))
