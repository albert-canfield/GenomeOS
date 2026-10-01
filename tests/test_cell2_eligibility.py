# SPDX-License-Identifier: AGPL-3.0-or-later
"""The independent-locus convention of genomeos/attribution/cell2.py, and its blinding.

Every test here pins a decision that was taken before any count: the 1 Mb span, the two joining
rules, the chaining, the conservative pooled count and the pre-set floor of 20. A run that loosened
any of them would fail these, which is the point of registering the convention in code.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from genomeos.attribution import cell2


@dataclass
class FakePair:
    """A pair as the counter is allowed to see it, plus fields it must never read."""

    cell: str
    chrom: str
    start: int
    end: int
    gene: str
    regulated: bool
    # fields outside BLINDED_FIELDS: a deletion value, a model score, anything predicted
    features: dict = field(default_factory=dict)
    prediction: float | None = None
    deletion: float | None = None


def key(chrom: str, start: int, gene: str, cell: str = "HCT116", length: int = 500) -> cell2.LocusKey:
    return cell2.LocusKey(cell, chrom, start, start + length, gene)


def test_span_and_floor_are_the_registered_figures():
    assert cell2.INDEPENDENT_LOCUS_SPAN == 1_000_000
    assert cell2.POOLED_LOCUS_FLOOR == 20
    assert cell2.PRIMARY_CELL == "K562"


def test_rule_text_states_both_joins_the_chaining_and_that_it_is_operational():
    rule = cell2.INDEPENDENT_LOCUS_RULE
    assert "share the measured gene" in rule
    assert "within 1 Mb" in rule
    assert "chained within a" in rule
    assert "not established biological independence" in rule


def test_blinded_fields_are_exactly_the_six_label_and_coordinate_fields():
    assert cell2.BLINDED_FIELDS == ("cell", "chrom", "start", "end", "gene", "regulated")
    assert set(cell2.LocusKey._fields) == {"cell", "chrom", "start", "end", "gene"}


def test_same_gene_joins_however_far_apart_the_elements_are():
    keys = [key("chr1", 1_000, "MYC"), key("chr1", 50_000_000, "MYC")]
    assert cell2.count_loci(keys) == 1


def test_elements_inside_the_span_join_although_the_genes_differ():
    keys = [key("chr2", 1_000_000, "A"), key("chr2", 1_900_000, "B")]
    assert cell2.gap(keys[0], keys[1]) == 899_500
    assert cell2.count_loci(keys) == 1


def test_elements_beyond_the_span_stay_separate():
    keys = [key("chr2", 1_000_000, "A"), key("chr2", 3_000_000, "B")]
    assert cell2.gap(keys[0], keys[1]) > cell2.INDEPENDENT_LOCUS_SPAN
    assert cell2.count_loci(keys) == 2


def test_the_span_boundary_is_inclusive_and_one_base_beyond_it_is_not():
    a = cell2.LocusKey("HCT116", "chr3", 0, 100, "A")
    at = cell2.LocusKey("HCT116", "chr3", 100 + cell2.INDEPENDENT_LOCUS_SPAN, 200 + 1_000_000, "B")
    beyond = cell2.LocusKey("HCT116", "chr3", 101 + cell2.INDEPENDENT_LOCUS_SPAN, 200 + 1_000_000, "B2")
    assert cell2.count_loci([a, at]) == 1
    assert cell2.count_loci([a, beyond]) == 2


def test_proximity_chains_within_a_chromosome_beyond_the_span():
    # each 0.9 Mb from the next, so the ends are 1.8 Mb apart and still one locus
    keys = [key("chr4", 0, "A"), key("chr4", 900_000, "B"), key("chr4", 1_800_000, "C")]
    assert cell2.count_loci(keys) == 1


def test_a_gap_breaks_the_chain_and_leaves_two_loci():
    keys = [
        key("chr4", 0, "A"),
        key("chr4", 900_000, "B"),
        key("chr4", 5_000_000, "C"),
        key("chr4", 5_900_000, "D"),
    ]
    assert cell2.count_loci(keys) == 2


def test_proximity_never_joins_across_chromosomes():
    keys = [key("chr5", 1_000, "A"), key("chr6", 1_000, "B")]
    assert cell2.count_loci(keys) == 2


def test_the_count_does_not_depend_on_the_input_order():
    keys = [
        key("chr7", 5_000_000, "C"),
        key("chr7", 0, "A"),
        key("chr8", 900_000, "B"),
        key("chr7", 900_000, "A2"),
    ]
    # chr7: 0 and 900,000 chain; 5,000,000 is 4.1 Mb from the chain's edge, so it stands alone. chr8
    # holds the third locus. Three either way round.
    assert cell2.count_loci(keys) == 3
    assert cell2.count_loci(list(reversed(keys))) == 3


def test_group_ids_partition_the_keys():
    keys = [key("chr9", 0, "A"), key("chr9", 900_000, "B"), key("chr9", 9_000_000, "C")]
    ids = cell2.group(keys)
    assert len(ids) == 3
    assert ids[0] == ids[1] != ids[2]


def test_overlapping_elements_have_no_gap():
    a = cell2.LocusKey("HCT116", "chr10", 1_000, 3_000, "A")
    b = cell2.LocusKey("HCT116", "chr10", 2_000, 4_000, "B")
    assert cell2.gap(a, b) == 0


# --- the shapes reported beside the count (reviewer's instruction, 2026-10-01) -----------------------


def test_locus_shapes_report_the_span_of_a_chained_locus_not_the_joining_span():
    keys = [key("chr4", 0, "A"), key("chr4", 900_000, "B"), key("chr4", 1_800_000, "C")]
    shapes = cell2.locus_shapes(keys)
    assert shapes["loci"] == 1
    assert shapes["span_mb"]["max"] == 1.8005  # 1,800,500 bp, wider than the 1 Mb join
    assert shapes["positives_per_locus"] == {"min": 3.0, "median": 3.0, "max": 3.0}


def test_locus_shapes_name_the_largest_locus_and_its_share():
    keys = [
        key("chr1", 0, "A"),
        key("chr1", 900_000, "B"),
        key("chr1", 1_800_000, "C"),
        key("chr2", 0, "D"),
    ]
    shapes = cell2.locus_shapes(keys)
    assert shapes["loci"] == 2
    assert shapes["largest_locus_positives"] == 3
    assert shapes["largest_locus_share_of_positives"] == 0.75
    assert shapes["largest_locus_span_mb"] == 1.8005
    assert shapes["positives_per_locus_counts"] == {"1": 1, "3": 1}
    assert shapes["span_mb"]["min"] == 0.0005


def test_locus_shapes_of_nothing_report_nulls_and_no_number():
    shapes = cell2.locus_shapes([])
    assert shapes["loci"] == 0
    assert shapes["span_mb"] == {"min": None, "median": None, "max": None}
    assert shapes["largest_locus_share_of_positives"] is None


def test_the_median_of_an_even_number_of_loci_is_the_mean_of_the_middle_two():
    keys = [key("chr1", 0, "A"), key("chr2", 0, "B", length=1_000_000)]
    shapes = cell2.locus_shapes(keys)
    assert shapes["span_mb"]["median"] == (0.0005 + 1.0) / 2


def heldout_like() -> list[FakePair]:
    """Two candidate cell types measuring the same two places, plus a K562 pair and a negative."""
    return [
        FakePair("HCT116", "chr1", 1_000, 1_500, "G1", True, prediction=0.9, deletion=-0.5),
        FakePair("HCT116", "chr1", 20_000_000, 20_000_500, "G2", True, prediction=0.1),
        FakePair("HCT116", "chr1", 30_000_000, 30_000_500, "G3", False, prediction=0.8),
        FakePair("Jurkat", "chr1", 1_000, 1_500, "G1", True, prediction=0.2),
        FakePair("K562", "chr1", 40_000_000, 40_000_500, "G4", True, prediction=0.7),
    ]


def test_eligibility_counts_per_cell_type_and_marks_the_candidates():
    out = cell2.eligibility(heldout_like(), ("K562", "HepG2", "GM12878", "IMR-90"))
    assert out["candidate_second_cell_types"] == ["HCT116", "Jurkat"]
    assert out["per_cell_type"]["HCT116"]["pairs"] == 3
    assert out["per_cell_type"]["HCT116"]["positives"] == 2
    assert out["per_cell_type"]["HCT116"]["independent_loci"] == 2
    assert out["per_cell_type"]["Jurkat"]["independent_loci"] == 1
    assert out["per_cell_type"]["K562"]["candidate_second_cell_type"] is False


def test_eligibility_carries_the_shapes_per_cell_type_and_pooled():
    out = cell2.eligibility(heldout_like(), ("K562", "HepG2", "GM12878", "IMR-90"))
    assert out["per_cell_type"]["HCT116"]["locus_shapes"]["loci"] == 2
    assert out["per_cell_type"]["HCT116"]["locus_shapes"]["largest_locus_positives"] == 1
    pooled = out["pooled_over_candidates"]["locus_shapes_genome_wide"]
    assert pooled["loci"] == 2
    # the chr1:1,000 place was measured by both candidate cell types: two positives, one locus
    assert pooled["largest_locus_positives"] == 2
    assert pooled["largest_locus_share_of_positives"] == round(2 / 3, 4)


def test_the_pooled_count_read_against_the_floor_is_the_conservative_one():
    out = cell2.eligibility(heldout_like(), ("K562", "HepG2", "GM12878", "IMR-90"))
    pooled = out["pooled_over_candidates"]
    assert pooled["positives"] == 3
    # HCT116 and Jurkat measured the same chr1:1,000 place: one locus genome-wide, two when summed
    assert pooled["independent_loci_genome_wide"] == 2
    assert pooled["independent_loci_summed_per_cell_type"] == 3
    assert pooled["which_is_read_against_the_floor"] == "independent_loci_genome_wide"
    assert out["floor"] == 20


def test_the_subset_with_a_deletion_value_is_counted_and_is_not_read_against_the_floor():
    pairs = heldout_like() + [
        FakePair("GM12878", "chr1", 1_000, 1_500, "G1", True),
        FakePair("GM12878", "chr5", 9_000_000, 9_000_500, "G9", True),
    ]
    out = cell2.eligibility(pairs, ("K562", "HepG2", "GM12878", "IMR-90"))
    subset = out["pooled_over_candidates_with_a_deletion_value"]
    assert subset["cell_types"] == ["GM12878"]
    assert subset["positives"] == 2
    assert subset["independent_loci_genome_wide"] == 2
    assert subset["read_against_the_floor"] is False
    # the floor verdict still reads the count over every candidate, which is the registered rule
    assert out["pooled_over_candidates"]["positives"] == 5
    assert (
        cell2.verdict(out)["pooled_independent_loci_genome_wide"]
        == (out["pooled_over_candidates"]["independent_loci_genome_wide"])
    )


def test_the_subset_is_empty_when_no_candidate_has_a_deletion_value():
    out = cell2.eligibility(heldout_like(), ("K562", "HepG2", "IMR-90"))
    subset = out["pooled_over_candidates_with_a_deletion_value"]
    assert subset["cell_types"] == []
    assert subset["positives"] == 0
    assert subset["independent_loci_genome_wide"] == 0
    assert "none" in subset["population"]


def test_deletion_availability_follows_the_model_cells_and_nothing_else():
    model = ("K562", "HepG2", "GM12878", "IMR-90")
    out = cell2.eligibility(heldout_like(), model)
    assert out["per_cell_type"]["K562"]["deletion_value_available"] is True
    assert out["per_cell_type"]["HCT116"]["deletion_value_available"] is False
    assert out["per_cell_type"]["Jurkat"]["deletion_value_available"] is False
    assert cell2.deletion_available("IMR-90", model) is True


def test_deletion_availability_matches_the_crispri_module_itself():
    from genomeos.attribution import crispri

    for cell in crispri.MODEL_CELLS:
        assert cell2.deletion_available(cell, crispri.MODEL_CELLS)
    for cell in ("HCT116", "WTC11", "Jurkat"):
        assert not cell2.deletion_available(cell, crispri.MODEL_CELLS)


def test_the_count_is_blind_to_every_field_outside_the_blinded_six():
    plain = heldout_like()
    loud = heldout_like()
    for i, p in enumerate(loud):
        p.prediction = float(i)
        p.deletion = -float(i)
        p.features = {"deletion_drop": -0.42 * i, "score": 7.0}
    model = ("K562", "HepG2", "GM12878", "IMR-90")
    assert cell2.eligibility(loud, model) == cell2.eligibility(plain, model)


def test_locus_key_carries_no_field_outside_the_blinded_six():
    p = heldout_like()[0]
    k = cell2.locus_key(p)
    assert k == cell2.LocusKey("HCT116", "chr1", 1_000, 1_500, "G1")
    assert 0.9 not in tuple(k)
    assert -0.5 not in tuple(k)


def test_the_verdict_applies_the_preset_floor_in_both_directions():
    counted = {"pooled_over_candidates": {"independent_loci_genome_wide": 19}}
    v = cell2.verdict(counted)
    assert v["floor"] == 20 and v["meets_floor"] is False and v["decision"] == "no-go, stop"
    at = cell2.verdict({"pooled_over_candidates": {"independent_loci_genome_wide": 20}})
    assert at["meets_floor"] is True and at["decision"] == "eligible to draft a registration"


def test_the_verdict_names_the_population_and_says_the_primary_cell_type_is_out_of_it():
    out = cell2.eligibility(heldout_like(), ("K562", "HepG2", "GM12878", "IMR-90"))
    v = cell2.verdict(out)
    assert v["cell_types_pooled"] == ["HCT116", "Jurkat"]
    assert "K562" not in v["cell_types_pooled"]
    assert v["primary_cell_type_excluded"] == "K562"
    assert v["primary_cell_type_is_in_the_pool"] is False
    assert "other than K562" in v["population"]
    assert v["positives_pooled"] == out["pooled_over_candidates"]["positives"]


def test_the_verdict_reports_what_carries_the_count():
    out = cell2.eligibility(heldout_like(), ("K562", "HepG2", "GM12878", "IMR-90"))
    carries = cell2.verdict(out)["what_carries_the_count"]
    shapes = out["pooled_over_candidates"]["locus_shapes_genome_wide"]
    assert carries["largest_locus_positives"] == shapes["largest_locus_positives"]
    assert carries["largest_locus_share_of_positives"] == shapes["largest_locus_share_of_positives"]
    assert carries["largest_locus_span_mb"] == shapes["largest_locus_span_mb"]
    assert carries["span_mb"] == shapes["span_mb"]
    assert carries["positives_per_locus"] == shapes["positives_per_locus"]
    assert carries["positives_per_locus_counts"] == shapes["positives_per_locus_counts"]


def test_the_shapes_beside_the_verdict_do_not_change_the_floor_decision():
    out = cell2.eligibility(heldout_like(), ("K562", "HepG2", "GM12878", "IMR-90"))
    v = cell2.verdict(out)
    pooled = out["pooled_over_candidates"]["independent_loci_genome_wide"]
    assert v["meets_floor"] is (pooled >= 20)
    assert v["floor"] == cell2.POOLED_LOCUS_FLOOR


# --- prior exposure, split under the pooled convention (reviewer, 2026-10-02) ------------------------


def exposure_pairs() -> list[FakePair]:
    """One locus held by an exposed and an unexposed cell type, one held by the unexposed one alone,
    one held by the exposed one alone, and a K562 pair that is none of the pool's business."""
    return [
        FakePair("GM12878", "chr1", 1_000, 1_500, "G1", True),
        FakePair("WTC11", "chr1", 1_200, 1_700, "G1b", True),
        FakePair("WTC11", "chr2", 5_000_000, 5_000_500, "G2", True),
        FakePair("HCT116", "chr3", 7_000_000, 7_000_500, "G3", True),
        FakePair("K562", "chr4", 1_000, 1_500, "G4", True),
        FakePair("WTC11", "chr2", 9_000_000, 9_000_500, "G5", False),
    ]


def test_the_exposed_cell_types_are_named_with_where_their_reading_sits():
    assert set(cell2.EXPOSED_CELL_TYPES) == {"GM12878", "HCT116"}
    for note in cell2.EXPOSED_CELL_TYPES.values():
        assert "crispri_published.json" in note


def test_exposure_splits_the_pooled_loci_and_excludes_the_primary_cell_type():
    e = cell2.exposure(exposure_pairs())
    assert e["pooled_positives"] == 4  # the K562 pair and the negative are out
    assert e["pooled_loci"] == 3
    assert e["already_exposed_loci"] == 2
    assert e["never_exposed_loci"] == 1
    assert e["never_exposed_cell_types"] == ["WTC11"]


def test_a_locus_shared_with_an_exposed_cell_type_is_counted_as_exposed():
    e = cell2.exposure(exposure_pairs())
    assert e["loci_in_more_than_one_cell_type"] == 1
    assert e["loci_counted_as_exposed_only_because_they_merge_with_an_exposed_cell_type"] == 1
    assert "counted as exposed" in e["merged_locus_side"]
    # the shared locus holds two positives and both are on the exposed side
    assert e["already_exposed_positives"] == 3
    assert e["never_exposed_positives"] == 1


def test_the_split_is_not_the_sum_of_the_per_cell_type_counts():
    e = cell2.exposure(exposure_pairs())
    r = e["reconciliation_of_the_summed_and_pooled_counts"]
    assert r["independent_loci_summed_per_cell_type"] == 4
    assert r["independent_loci_pooled_genome_wide"] == 3
    assert r["difference"] == 1
    assert r["pooled_loci_that_absorb_more_than_one_per_cell_type_locus"] == 1
    assert e["already_exposed_loci"] + e["never_exposed_loci"] == e["pooled_loci"]


def test_the_never_scored_stratum_is_grouped_on_its_own_and_marked_descriptive():
    own = cell2.exposure(exposure_pairs())["never_scored_stratum_on_its_own"]
    assert own["positives"] == 2  # both WTC11 positives, including the one sharing an exposed locus
    assert own["independent_loci_genome_wide"] == 2
    assert "descriptive only" in own["status"]
    assert "WTC11" in own["population"]


def test_exposure_takes_an_explicit_exposed_set_without_touching_the_registered_one():
    e = cell2.exposure(exposure_pairs(), exposed={"WTC11": "for the test only"})
    assert e["never_exposed_cell_types"] == ["GM12878", "HCT116"] or set(e["never_exposed_cell_types"]) <= {
        "GM12878",
        "HCT116",
    }
    assert set(cell2.EXPOSED_CELL_TYPES) == {"GM12878", "HCT116"}


# --- the power gate's registered line, readings and profile feasibility ------------------------------


def test_the_power_gate_line_and_design_are_registered_figures():
    assert cell2.POWER_GATE_SHARE_FLOOR == 0.5
    assert "share of subsamples whose interval excludes zero" in cell2.POWER_GATE_DESIGN
    assert "positives-per-locus profile" in cell2.POWER_GATE_DESIGN
    assert "not to spend" in cell2.POWER_GATE_DESIGN


def test_the_readings_are_fixed_in_advance_and_refuse_no_effect():
    assert cell2.READINGS["lower_bound_above_zero"] == "replicated outside K562 on these cell types"
    assert cell2.READINGS["interval_covering_zero"] == "no difference detected"
    assert "no effect" in cell2.READINGS["never_say"]


def test_a_profile_is_matchable_only_by_loci_of_exactly_the_required_sizes():
    source = (
        [key("chr1", i * 10_000_000, f"A{i}") for i in range(3)]
        + [key("chr2", 0, "B1"), key("chr2", 100, "B2")]
        + [key("chr3", 0, "C1"), key("chr3", 100, "C2"), key("chr3", 200, "C3")]
    )
    # source holds three loci of one positive, one of two and one of three
    assert cell2._count_sizes(source) == {1: 3, 2: 1, 3: 1}
    ok = cell2.profile_match({"1": 2, "3": 1}, source)
    assert ok["exact_profile_available"] is True
    assert ok["shortfall_by_locus_size"] == {}


def test_a_profile_needing_a_larger_locus_than_the_source_holds_is_refused_with_the_shortfall():
    source = [key("chr1", 0, "A1"), key("chr1", 100, "A2")]  # one locus of two positives
    out = cell2.profile_match({"2": 1, "5": 1}, source)
    assert out["exact_profile_available"] is False
    assert out["shortfall_by_locus_size"] == {"5": {"required": 1, "available": 0, "missing": 1}}
    assert "must be registered before it is used" in out["rule"]


def test_a_larger_source_locus_cannot_stand_in_for_a_smaller_required_one():
    source = [key("chr1", 0, "A1"), key("chr1", 100, "A2"), key("chr1", 200, "A3")]  # one locus of 3
    out = cell2.profile_match({"1": 1}, source)
    assert out["exact_profile_available"] is False
    assert out["shortfall_by_locus_size"]["1"]["missing"] == 1
