# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-pbrun: the counted clause, the carried minimum, and the things that must not drift.

These tests hold the four claims the lane's value rests on. (1) The minimum is lane-pseudobulk's
OWN dict, imported and not copied, so the two lanes cannot drift apart on either number. (2) M2h's
scale enters only as the fixed constant of `ac5580b` and NOTHING in this lane fits anything -- the
module's own syntax tree is read and a fitter in it fails this file. (3) Neither clause of the
minimum rescues the other: four cell types with 999 genes is not a qualifying inventory and thirteen
cell types would not be either. (4) The one deviation from M2h -- centring over the four compared
cell types rather than thirteen -- is applied to BOTH arms, and a planted test shows what centring
them differently would do.
"""

from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest
import tracked_paths

from genomeos.attribution import pbrun as pr
from genomeos.attribution import pseudobulk as pb

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "genomeos" / "attribution" / "pbrun.py"

#: A name whose presence in this lane's module would mean something is being fitted here.
FITTERS = frozenset(
    {
        "polyfit",
        "lstsq",
        "curve_fit",
        "minimize",
        "least_squares",
        "LinearRegression",
        "fit",
        "fit_transform",
        "optimize",
        "nnls",
    }
)


def test_the_minimum_is_lane_pseudobulks_own_dict_and_not_a_copy() -> None:
    assert pr.MINIMUM is pb.MINIMUM
    assert pr.MINIMUM["cell_types"] == 4
    assert pr.MINIMUM["genes"] == 1000
    assert pb.MINIMUM["genes_clause"] == pr.GENE_CLAUSE
    assert pr.WORDING_LIMIT == pb.WORDING_LIMIT


def test_the_scale_is_the_fixed_constant_of_ac5580b_and_m2hs_committed_value() -> None:
    tracked_paths.must_be_present(
        pr.MAGCAL_RESULT, was="M2h's committed record, which fixes the constant this lane may use"
    )
    assert pr.SCALE is pb.SCALE
    assert pr.SCALE == 0.2145911381525326
    committed = json.loads(pr.MAGCAL_RESULT.read_text())["per_form"]["multiplicative"]["scale"]
    assert committed == pr.SCALE


def test_no_fitter_appears_anywhere_in_this_lanes_module() -> None:
    tree = ast.parse(MODULE.read_text())
    used: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute):
            used.add(node.attr)
        elif isinstance(node, ast.Name):
            used.add(node.id)
    assert not (used & FITTERS), f"a fitter appears in {MODULE.name}: {sorted(used & FITTERS)}"


def test_the_scale_is_assigned_from_lane_pseudobulk_and_from_nothing_else() -> None:
    tree = ast.parse(MODULE.read_text())
    assignments = [
        node
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "SCALE" for t in node.targets)
    ]
    assert len(assignments) == 1
    value = assignments[0].value
    assert isinstance(value, ast.Attribute)
    assert value.attr == "SCALE"
    assert isinstance(value.value, ast.Name)
    assert value.value.id == "pb"


def test_neither_clause_of_the_minimum_rescues_the_other() -> None:
    assert pr.verdict(4, 1000) == "QUALIFIES"
    assert pr.verdict(4, 999) == "NO_SOURCE_QUALIFIES"
    assert pr.verdict(13, 999) == "NO_SOURCE_QUALIFIES"
    assert pr.verdict(3, 10_000) == "NO_SOURCE_QUALIFIES"
    assert pr.verdict(4, None) == "NO_SOURCE_QUALIFIES"
    assert pr.verdict(4, 1000) == pb.verdict(4, 1000)


def test_the_four_is_arithmetic_about_the_registered_control_and_not_prose() -> None:
    alpha = 0.05
    attainable = [n for n in range(1, 8) if math.factorial(n) >= 1 / alpha]
    assert min(attainable) == pr.MINIMUM["cell_types"]
    assert math.factorial(3) < 1 / alpha <= math.factorial(4)


def test_a_par_y_row_is_never_equal_to_the_plain_gene_it_shadows() -> None:
    assert pr.unversioned("ENSG00000182378.14") == "ENSG00000182378"
    assert pr.unversioned("ENSG00000182378.14_PAR_Y") == "ENSG00000182378.14_PAR_Y"
    assert pr.unversioned("ENSG00000182378.14_PAR_Y") != pr.unversioned("ENSG00000182378.14")


def test_resolves_demands_exactly_one_row_in_every_contributing_file() -> None:
    one = {"ENSG1": 1, "ENSG2": 1}
    duplicated = {"ENSG1": 2, "ENSG2": 1}
    missing = {"ENSG2": 1}
    assert pr.resolves([one, one], "ENSG1")
    assert not pr.resolves([one, duplicated], "ENSG1")
    assert not pr.resolves([one, missing], "ENSG1")
    assert not pr.resolves([], "ENSG1")


def test_the_held_out_population_reproduces_m2hs_committed_gene_count() -> None:
    tracked_paths.must_be_present(
        pr.MULTRNA_RESULT, was="7h's committed genes table, the population M2h calibrated"
    )
    tracked_paths.must_be_present(pr.MAGCAL_RESULT, was="M2h's committed record")
    population = pr.held_out_population()
    committed = json.loads(pr.MAGCAL_RESULT.read_text())["per_form"]["multiplicative"]
    assert len(population) == committed["split"]["held_out_genes"]
    assert len({g["symbol"] for g in population}) == len(population)


def test_the_split_is_m2hs_and_puts_a_gene_wholly_on_one_side() -> None:
    expected = int(hashlib.sha256(b"TP53").hexdigest()[:8], 16) % 2 == 1
    assert pr.is_held_out("TP53") == expected
    assert pr.is_held_out("TP53") == pr.is_held_out("TP53")


def test_centring_is_over_the_compared_cell_types_and_the_same_for_both_arms() -> None:
    four = np.array([[1.0, 2.0, 3.0, 4.0], [0.0, 0.0, 0.0, 4.0]])
    d = pr.centred(four)
    assert np.allclose(d.mean(axis=1), 0.0)
    thirteen = np.concatenate([four, np.zeros((2, 9))], axis=1)
    over_thirteen = pr.centred(thirteen)[:, :4]
    assert not np.allclose(d, over_thirteen), (
        "centring one arm over thirteen cell types and the other over four would make the two "
        "numbers incomparable; this is the error pr.centred exists to prevent"
    )


def test_delta_mse_is_m2hs_closed_form_at_a_fixed_scale_and_refits_nothing() -> None:
    rng = np.random.default_rng(7)
    d = rng.normal(size=(50, 4))
    p = rng.normal(size=(50, 4))
    a = pr.SCALE
    expected = float(np.mean(d**2) - np.mean((d - a * p) ** 2))
    assert pr.delta_mse(d, p, a) == pytest.approx(expected, abs=1e-12)
    best = float(np.sum(d * p) / np.sum(p**2))
    assert pr.delta_mse(d, p, best) >= pr.delta_mse(d, p, a)
    assert pr.delta_mse(d, p, a) != pr.delta_mse(d, p, best)


def test_the_pseudobulk_cpm_divides_by_the_sum_over_every_row_not_the_counted_ones() -> None:
    counts = np.array([1.0, 3.0, 996.0])
    cpm = pr.pseudobulk_cpm(counts)
    assert cpm.sum() == pytest.approx(1e6)
    assert cpm[0] == pytest.approx(1000.0)
    assert pr.pseudobulk_cpm(counts)[1] != pr.pseudobulk_cpm(counts[:2])[1]
    with pytest.raises(ValueError):
        pr.pseudobulk_cpm(np.zeros(3))


def test_the_measured_transform_carries_m2hs_pseudocount_of_one() -> None:
    assert pr.log2p1(np.array([0.0, 1.0, 3.0])).tolist() == [0.0, 1.0, 2.0]


def test_the_rsem_reader_counts_rows_per_unversioned_gene_id(tmp_path: Path) -> None:
    table = tmp_path / "t.tsv"
    table.write_text(
        "gene_id\tlength\texpected_count\tTPM\n"
        "ENSG1.1\t100\t5\t1.0\n"
        "ENSG2.3\t100\t0\t0.0\n"
        "ENSG2.4\t100\t0\t0.0\n"
        "ENSG3.1_PAR_Y\t100\t0\t0.0\n"
        "gSpikein1\t100\t0\t0.0\n"
    )
    counts = pr.row_counts_rsem_tsv(table)
    assert counts["ENSG1"] == 1
    assert counts["ENSG2"] == 2
    assert "ENSG3" not in counts
    assert counts["ENSG3.1_PAR_Y"] == 1


def test_an_archive_with_two_feature_tables_raises_rather_than_guessing() -> None:
    with pytest.raises(ValueError):
        pr._one_member(["a/features.tsv", "b/features.tsv"], "features.tsv")
    with pytest.raises(ValueError):
        pr._one_member(["a/barcodes.tsv"], "features.tsv")
    assert pr._one_member(["a/features.tsv", "a/barcodes.tsv"], "features.tsv") == "a/features.tsv"


def test_the_four_cell_types_are_read_from_lane_scsweeps_record_and_not_hardcoded() -> None:
    tracked_paths.must_be_present(
        pr.SCSWEEP_RESULT, was="lane-scsweep's committed record, which names the four cell types"
    )
    cells = pr.qualifying_cells()
    assert len(cells) == pr.MINIMUM["cell_types"]
    record = json.loads(pr.SCSWEEP_RESULT.read_text())
    assert record["exact_label_count"] == len(cells)
    assert "CD14_positive_monocyte" not in cells, (
        "the fifth source clears every clause but the label string, so it is in "
        "same_biosample_count and NOT in this lane's four"
    )


def test_the_aggregation_states_the_sum_the_normalisation_and_the_absence_of_a_filter() -> None:
    assert "SUM" in pr.AGGREGATION["statement"]
    assert "counts per million" in pr.AGGREGATION["statement"]
    assert "PSEUDOCOUNT 1.0" in pr.AGGREGATION["statement"]
    assert "NO FILTER AT ANY STEP" in pr.AGGREGATION["no_filter"]
    assert "cpm_not_tpm" in pr.AGGREGATION


def test_the_modality_of_the_four_is_declared_before_any_value_is_read() -> None:
    assert "ONE IS SINGLE-NUCLEUS" in pr.SC_SELECTION["modality_note"]
    assert "ovary" in pr.SC_SELECTION["modality_note"]


def test_the_form_is_the_one_m2h_read_better_on() -> None:
    tracked_paths.must_be_present(pr.MAGCAL_RESULT, was="M2h's committed record")
    per_form = json.loads(pr.MAGCAL_RESULT.read_text())["per_form"]
    better = [name for name, f in per_form.items() if f["outcome"] == "BETTER"]
    assert better == [pr.FORM]


# --- the committed count record, and the registration this lane may now write --------------------
#
# These read `data/results/pbrun_gene_count.json`, which `3ef7e82` makes git track, so they assert
# it with `tracked_paths.must_be_present` and NEVER with a skip: a skip keyed on a tracked path can
# never fire and reads to a later reviewer as a condition that is a constant.


def test_the_committed_count_clears_both_numbers_of_the_minimum() -> None:
    tracked_paths.must_be_present(
        pr.GENE_COUNT_PATH, was="this lane's own committed count of the 1,000-gene clause"
    )
    count = json.loads(pr.GENE_COUNT_PATH.read_text())
    assert count["qualifying_cell_count"] >= pr.MINIMUM["cell_types"]
    assert count["qualifying_gene_count"] >= pr.MINIMUM["genes"]
    assert count["answer"] == pr.verdict(count["qualifying_cell_count"], count["qualifying_gene_count"])
    assert count["answer"] == "QUALIFIES"


def test_the_counted_population_is_a_subset_of_m2hs_held_out_side() -> None:
    tracked_paths.must_be_present(pr.GENE_COUNT_PATH, was="this lane's committed count")
    count = json.loads(pr.GENE_COUNT_PATH.read_text())
    genes = count["qualifying_genes"]
    assert len(genes) == len(set(genes)) == count["qualifying_gene_count"]
    assert all(pr.is_held_out(symbol) for symbol in genes)
    assert len(genes) <= count["held_out_population"]


def test_no_cell_type_is_counted_on_fewer_genes_than_the_shared_population() -> None:
    tracked_paths.must_be_present(pr.GENE_COUNT_PATH, was="this lane's committed count")
    count = json.loads(pr.GENE_COUNT_PATH.read_text())
    per_cell = {e["cell"]: e["usable_on_both_sides"] for e in count["per_cell"]}
    assert set(per_cell) == set(count["qualifying_cells"])
    assert min(per_cell.values()) >= count["qualifying_gene_count"], (
        "the shared population is the INTERSECTION over the four cell types, so it cannot exceed "
        "any one of them"
    )


def test_the_count_record_names_every_input_it_read_including_the_two_h5ad_files() -> None:
    tracked_paths.must_be_present(pr.GENE_COUNT_PATH, was="this lane's committed count")
    manifest = json.loads(pr.GENE_COUNT_PATH.read_text())["result_manifest"]
    paths = {i["path"] for i in manifest["inputs"]}
    assert manifest["input_count"] == len(manifest["inputs"])
    assert sum(1 for p in paths if p.endswith(".h5ad")) == 2, (
        "both h5ad matrices are opened with h5py, which the read tracer does not see; they are "
        "declared by hand and a missing one makes the result unrebuildable (9ce4c6e)"
    )
    assert all(not Path(p).is_absolute() for p in paths), (
        "input_entry records the path verbatim, so an absolute path makes the result unrebuildable (d72f959)"
    )
    assert manifest["code_cleanliness"]["own_code_is_committed"]


def test_the_margin_is_half_of_m2hs_committed_held_out_gain_and_not_a_round_number() -> None:
    tracked_paths.must_be_present(pr.MAGCAL_RESULT, was="M2h's committed record")
    committed = json.loads(pr.MAGCAL_RESULT.read_text())["per_form"][pr.FORM]
    assert pr.margin() == committed["held_out_delta_mse"] / 2.0
    assert pr.margin() > 0.0
    assert "half of M2h's committed held-out Delta MSE" in pr.MARGIN_RULE


def test_the_statistic_is_m2hs_unchanged() -> None:
    assert pr.BOOTSTRAP["resamples"] == 2000
    assert pr.BOOTSTRAP["seed"] == 20261010
    assert pr.PERMUTATIONS == {"count": 1000, "seed": 20261011}
    assert pr.CARRIED["pseudocount"] == 1.0
    assert pr.CARRIED["min_genes"] == 100


def test_the_three_outcomes_are_declared_and_inconclusive_licenses_nothing() -> None:
    assert set(pr.DECISION) == {
        "agreement",
        "disagreement",
        "inconclusive",
        "falsifier",
        "stop_rule",
    }
    assert "FOR THIS ENDPOINT ONLY" in pr.DECISION["agreement"]
    assert "NEVER 'single-nucleus RNA equals bulk RNA'" in pr.DECISION["agreement"]
    assert "CLOSES WITHOUT A LIFT" in pr.DECISION["disagreement"]
    assert "licenses nothing" in pr.DECISION["inconclusive"]
    assert "weak agreement" in pr.DECISION["inconclusive"]


def test_the_exposure_of_the_bulk_side_is_stated_and_not_softened() -> None:
    assert "already read by items 7, 7d, 7h, M2h and M2i" in pr.EXPOSURE
    assert "NOT A FRESH MEASUREMENT" in pr.EXPOSURE


def test_the_one_deviation_from_m2h_is_named_as_a_deviation() -> None:
    assert "ONE, NAMED AS A DEVIATION" in pr.DEVIATION_FROM_M2H
    assert "THIRTEEN" in pr.DEVIATION_FROM_M2H
    assert "FOUR" in pr.DEVIATION_FROM_M2H


def test_the_limitations_name_the_single_nucleus_count_and_the_normalisation() -> None:
    joined = " ".join(pr.LIMITATIONS)
    assert "ONLY ONE OF THE FOUR CELL TYPES IS SINGLE-NUCLEUS" in joined
    assert "CPM" in joined
    assert "length-normalised TPM" in joined
    assert "FOUR CELL TYPES, not thirteen" in joined


# --- the decision, as code rather than as prose in a script --------------------------------------


def test_an_arms_reading_is_m2hs_rule_with_the_straddle_branches_collapsed() -> None:
    assert pr.arm_reads([0.01, 0.05]) == "better"
    assert pr.arm_reads([-0.05, -0.01]) == "worse"
    assert pr.arm_reads([-0.01, 0.05]) == "inconclusive"
    assert pr.arm_reads([0.0, 0.05]) == "inconclusive"
    assert pr.arm_reads([-0.05, 0.0]) == "inconclusive"


def test_agreement_needs_both_clauses_and_nothing_less() -> None:
    m = 0.02
    assert pr.decide([0.03, 0.05], [0.0, 0.01], m) == "AGREEMENT"
    assert pr.decide([0.03, 0.05], [-0.02, 0.02], m) == "AGREEMENT"
    # better, but the difference is not resolved against the margin
    assert pr.decide([0.03, 0.05], [0.01, 0.03], m) == "INCONCLUSIVE"
    # the difference is inside the margin, but the arm itself does not read better
    assert pr.decide([-0.01, 0.05], [0.0, 0.01], m) == "INCONCLUSIVE"


def test_the_falsifier_takes_precedence_over_everything_else() -> None:
    m = 0.02
    assert pr.decide([-0.05, -0.01], [0.0, 0.001], m) == "DISAGREEMENT"
    assert pr.decide([0.03, 0.05], [0.03, 0.06], m) == "DISAGREEMENT"
    assert pr.decide([0.03, 0.05], [-0.06, -0.03], m) == "DISAGREEMENT"


def test_the_decision_has_exactly_three_outcomes_over_a_grid() -> None:
    m = pr.margin()
    outcomes = set()
    for lo in (-0.05, -0.001, 0.0, 0.001, 0.05):
        for hi in (lo, lo + 0.001, lo + 0.1):
            for dlo in (-0.1, -m, 0.0, m, 0.1):
                for dhi in (dlo, dlo + 0.0005, dlo + 0.2):
                    outcomes.add(pr.decide([lo, hi], [dlo, dhi], m))
    assert outcomes == {"AGREEMENT", "DISAGREEMENT", "INCONCLUSIVE"}


def test_the_within_gene_permutation_keeps_every_genes_own_values() -> None:
    p = np.array([[1.0, 2.0, 3.0, 4.0], [10.0, 20.0, 30.0, 40.0]])
    draws = list(pr.permute_within_gene(p, seed=pr.PERMUTATIONS["seed"], count=25))
    assert len(draws) == 25
    for permuted in draws:
        assert np.array_equal(np.sort(permuted, axis=1), np.sort(p, axis=1))
    assert any(not np.array_equal(permuted, p) for permuted in draws), (
        "a control that never moves anything would be no control"
    )
    again = list(pr.permute_within_gene(p, seed=pr.PERMUTATIONS["seed"], count=25))
    assert all(np.array_equal(a, b) for a, b in zip(draws, again, strict=True))


def test_the_control_p_is_m2hs_formula() -> None:
    permuted = np.array([0.0, 0.1, 0.2, 0.3])
    assert pr.control_p(0.4, permuted) == 1 / 5
    assert pr.control_p(0.0, permuted) == 5 / 5


def test_the_interval_is_the_9833_percent_percentile_interval() -> None:
    values = np.arange(10_000, dtype=np.float64)
    low, high = pr.percentile_interval(values)
    assert low == pytest.approx(np.percentile(values, 0.8333))
    assert high == pytest.approx(np.percentile(values, 99.1667))
