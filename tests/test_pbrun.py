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
