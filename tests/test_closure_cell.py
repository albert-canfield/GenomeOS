# SPDX-License-Identifier: AGPL-3.0-or-later
"""The cell-level closure (M3): the sign convention, the arms, the layout and the decision rule.

Every number this test asserts is one a defect would move. The two that matter most are the SIGN
(an element that holds a gene up must give a positive effect, and a depleted element must give a
negative readout -- a sign slip here would turn the registered falsifier into its opposite) and the
DECISION RULE, which is exercised on each of its five branches with tables built by hand, so the
rule the registration fixed is the rule the run applies.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos.attribution import closure as gene_level  # noqa: E402
from genomeos.attribution import closure_cell as cc  # noqa: E402

GUIDE_BED_HOW = (
    "data/knowledge is a local store: fetch the screen library's design with "
    "scripts/closure_cell_inventory.py, which downloads ENCODE ENCFF630CYG"
)


def _script(name: str):
    """One of this lane's scripts, loaded as a module; they are scripts, as their siblings are."""
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# ---- the sign convention -------------------------------------------------------------------------


def test_the_cell_s_effect_is_signed_the_way_the_gene_level_closure_signs_it():
    """A fall in the gene on deleting the element means the element HOLDS THE GENE UP: positive."""
    element = {"predicted_coding_by_cell": {"K562": -0.42, "HepG2": 0.1}}
    assert cc.element_effect(element) == pytest.approx(0.42)
    assert cc.element_effect(element) == pytest.approx(gene_level.effect_of(-0.42))
    assert cc.element_effect({"predicted_coding_by_cell": {"K562": 0.3}}) == pytest.approx(-0.3)
    assert cc.element_effect({"predicted_coding_by_cell": {"HepG2": 0.3}}) is None
    assert cc.element_effect({}) is None


def test_a_depleted_element_reads_negative_and_an_enriched_one_positive():
    depleted = cc.log2_ratio(10.0, 1_000_000.0, 100.0, 1_000_000.0)
    enriched = cc.log2_ratio(100.0, 1_000_000.0, 10.0, 1_000_000.0)
    assert depleted is not None and depleted < 0
    assert enriched is not None and enriched > 0
    assert depleted == pytest.approx(-enriched)
    # an element the screen emptied keeps a finite value instead of leaving its arm
    assert cc.log2_ratio(0.0, 1_000_000.0, 500.0, 1_000_000.0) == pytest.approx(math.log2(1.0 / 501.0))
    assert cc.log2_ratio(5.0, 0.0, 5.0, 10.0) is None


# ---- overlap and the nearest coding gene ---------------------------------------------------------


def test_overlap_is_half_open_and_crosses_a_bin_boundary():
    rows = [(9_990, 10_010, "spanning"), (20_000, 20_003, "inside"), (30_000, 30_003, "far")]
    index = cc.bin_index(rows)
    assert [rows[i][2] for i in cc.overlapping(rows, index, 9_995, 9_998)] == ["spanning"]
    assert [rows[i][2] for i in cc.overlapping(rows, index, 10_005, 10_008)] == ["spanning"]
    # half-open on both sides: a guide that ends where the element starts does not overlap it
    assert cc.overlapping(rows, index, 19_997, 20_000) == []
    assert [rows[i][2] for i in cc.overlapping(rows, index, 20_002, 20_005)] == ["inside"]
    # an element present in two bins is returned once, not twice
    assert len(cc.overlapping(rows, index, 9_990, 10_010)) == 1


def test_the_nearest_coding_tss_is_zero_inside_the_element_and_the_closer_side_outside():
    positions = [1_000, 5_000, 9_000]
    symbols = ["LEFT", "MIDDLE", "RIGHT"]
    assert cc.nearest_coding(positions, symbols, 4_900, 5_100) == ("MIDDLE", 0)
    assert cc.nearest_coding(positions, symbols, 5_200, 5_300) == ("MIDDLE", 200)
    # the distance is counted from the element's nearest BASE, so 8_099 to 9_000 is 901
    assert cc.nearest_coding(positions, symbols, 8_000, 8_100) == ("RIGHT", 901)
    assert cc.nearest_coding([], [], 10, 20)[0] is None


# ---- the dependency call and the arms ------------------------------------------------------------


def test_the_dependency_call_has_a_band_that_is_in_no_arm():
    effects = {"HARD": -0.9, "AT_THE_BAR": cc.ESSENTIAL, "BORDER": -0.3, "FREE": 0.05}
    assert cc.dependency(effects, "HARD") == "essential"
    assert cc.dependency(effects, "AT_THE_BAR") == "essential"  # at or below the bar
    assert cc.dependency(effects, "BORDER") == "borderline"
    assert cc.dependency(effects, "FREE") == "not"
    assert cc.dependency(effects, "ABSENT") == "unmeasured"
    assert cc.dependency(effects, None) == "unmeasured"


def test_the_arms_separate_the_attribution_from_the_neighbourhood():
    assert cc.arm_of("essential", "not", 0.5) == "attributed_essential_neighbourhood_not"
    assert cc.arm_of("not", "not", 0.5) == "attributed_not_neighbourhood_not"
    assert cc.arm_of("not", "essential", 0.5) == "attributed_not_neighbourhood_essential"
    assert cc.arm_of("essential", "essential", 0.5) == "attributed_essential_neighbourhood_essential"
    assert set(cc.ARMS) == {cc.arm_of(a, b, 0.5) for a in ("essential", "not") for b in ("essential", "not")}
    # a prediction under the magnitude bar is in no arm, and the opposite sign goes to its own
    assert cc.arm_of("essential", "not", cc.MIN_EFFECT / 2) is None
    assert cc.arm_of("essential", "not", -0.5) == cc.DIRECTION_ARM
    assert cc.arm_of("not", "not", -0.5) is None
    assert cc.arm_of("borderline", "not", 0.5) is None
    assert cc.arm_of("essential", "unmeasured", 0.5) is None
    assert cc.DIRECTION_ARM not in cc.ARMS


# ---- the published layout ------------------------------------------------------------------------


#: One row as ENCODE publishes it, taken from the first line of ENCFF995GEZ: eighteen fields, a
#: 3 bp target interval, the replicate's count at field 5 and the guide type at field 16.
A_PUBLISHED_ROW = (
    "chr1\t605582\t605585\tNA\tchr1:605582-605585:+\t446\t+\tchr1:605582-605585:+\t"
    "NA\tNA\tNA\tNA\tNA\tNA\tAACCACAGGGAGAGCTGCCG\tGAACCACAGGGAGAGCTGCCG\ttargeting\tNA"
)


def test_the_published_row_parses_and_a_malformed_one_is_refused_rather_than_guessed_at():
    row = cc.parse_guide_row(A_PUBLISHED_ROW)
    assert row == ("chr1", 605582, 605585, 446.0, "targeting")
    fields = A_PUBLISHED_ROW.split("\t")
    assert len(fields) == cc.FIELD_COUNT
    assert fields[cc.FIELDS["count"]] == "446"
    assert fields[cc.FIELDS["guide_type"]] == cc.TARGETING
    assert cc.parse_guide_row("\t".join(fields[:-1])) is None  # a field short
    assert cc.parse_guide_row(A_PUBLISHED_ROW.replace("\t446\t", "\tNA\t")) is None
    assert cc.parse_guide_row(A_PUBLISHED_ROW.replace("chr1\t", "1\t", 1)) is None
    assert cc.parse_guide_row(A_PUBLISHED_ROW.replace("605582\t605585", "605585\t605582", 1)) is None


def test_every_row_is_counted_somewhere_and_the_totals_close():
    rows = [(605_500, 605_600, "EH_A")]
    index = {"chr1": (rows, cc.bin_index(rows))}
    per_element: dict[str, float] = {}
    guides: dict[str, int] = {}
    census = cc.new_census()
    non_targeting = {"guides": 0.0, "count": 0.0}
    control_row = A_PUBLISHED_ROW.replace("\ttargeting\t", "\tnegative_control\t").replace(
        "\t446\t", "\t50\t"
    )
    off = A_PUBLISHED_ROW.replace("chr1\t605582\t605585", "chr1\t900000\t900003", 1).replace(
        "\t446\t", "\t7\t"
    )
    unplaced = A_PUBLISHED_ROW.replace("chr1\t", "chrUn_KI270742v1\t", 1).replace("\t446\t", "\t3\t")
    for line in (A_PUBLISHED_ROW, control_row, off, unplaced, "\t".join(["x"] * 3), ""):
        cc.fold_row(line, index, per_element, guides, census, non_targeting)
    assert census["rows"] == 5  # the empty line is no row
    assert census["malformed"] == 1
    assert census["targeting"] == 3
    assert census["non_targeting"] == 1
    assert census["outside_the_population"] == 1
    assert census["unplaced_chromosome"] == 1
    assert per_element == {"EH_A": 446.0}
    assert guides == {"EH_A": 1}
    assert non_targeting == {"guides": 1.0, "count": 50.0}
    assert census["count_outside_the_population"] == 10.0
    total = sum(per_element.values()) + non_targeting["count"] + census["count_outside_the_population"]
    assert total == 446.0 + 50.0 + 7.0 + 3.0


# ---- the statistics ------------------------------------------------------------------------------


def test_the_rank_test_finds_a_shift_and_reports_nothing_when_an_arm_is_empty():
    below = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3]
    floor = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    separated = cc.mann_whitney_u(below, floor)
    assert separated["u"] == 0  # every value of the first arm ranks below every value of the second
    assert separated["p_one_sided_a_below_b"] is not None
    assert separated["p_one_sided_a_below_b"] < 0.001
    # the same two arms the other way round must NOT read as depletion
    assert cc.mann_whitney_u(floor, below)["p_one_sided_a_below_b"] > 0.999
    same = cc.mann_whitney_u(floor, list(floor))
    assert same["p_one_sided_a_below_b"] == pytest.approx(0.5, abs=0.05)
    assert cc.mann_whitney_u([], floor)["p_one_sided_a_below_b"] is None
    # all ties: the tie correction leaves no variance, so no p-value is invented
    assert cc.mann_whitney_u([1.0, 1.0], [1.0, 1.0])["p_one_sided_a_below_b"] is None


def test_the_depletion_cut_makes_the_floor_arm_s_rate_five_percent_by_construction():
    floor = [float(i) for i in range(100)]
    cut = cc.percentile(floor, cc.DEPLETION_PERCENTILE)
    assert cut is not None
    assert sum(1 for v in floor if v <= cut) == 5
    assert cc.percentile([], 5) is None
    assert cc.percentile([3.0], 5) == 3.0
    assert cc.median([1.0, 2.0, 3.0, 4.0]) == pytest.approx(2.5)
    assert cc.mean([]) is None


# ---- the registration's own guard ----------------------------------------------------------------


def test_the_registration_refuses_a_figure_the_run_would_fill_at_any_depth():
    register = _script("closure_cell_register")
    clean = {"question": "x", "arms": {"sizes_known_in_advance": {"a": 1}}}
    assert register.forbidden_at_any_depth(clean) == []
    nested = {"arms": {"table": [{"median_a": -0.2}]}}
    assert register.forbidden_at_any_depth(nested) == ["arms.table[0].median_a"]
    assert "outcome" in register.FORBIDDEN
    assert "p_ad" in register.FORBIDDEN
    for wording in register.OUTCOMES.values():
        assert "{" in wording


def test_the_registered_wordings_state_a_denominator_and_never_claim_equality():
    register = _script("closure_cell_register")
    for name, wording in register.OUTCOMES.items():
        assert "{n_" in wording, name
        assert " equal" not in wording.lower(), name
    assert "not no difference" in register.OUTCOMES["no_difference_detected"]


# ---- the decision rule, branch by branch ---------------------------------------------------------


def _table(**sizes_and_medians):
    """A minimal arm table: {arm: (n, median)}, shaped as the run's own `arm_table` shapes it."""
    return {
        arm: {
            "measured_elements": n,
            "median_log2_ratio": median,
            "rate_with_denominator": f"3 of {n}",
        }
        for arm, (n, median) in sizes_and_medians.items()
    }


def _tests(**p_by_arm):
    return {
        f"{arm}_below_attributed_not_neighbourhood_not": {"p_one_sided_a_below_b": p}
        for arm, p in p_by_arm.items()
    }


def test_the_decision_rule_selects_each_outcome_on_the_table_that_defines_it():
    run = _script("closure_cell_run")
    floor = {run.ARM_D: (200, 0.0)}
    working = {cc.POSITIVE_CONTROL: (50, -1.0)}

    stopped = run.select(_table(**{**floor, run.ARM_A: (5, -0.5)}, **working), _tests(), 30)
    assert stopped[0] == "stopped_for_power"

    broken = _table(**{**floor, run.ARM_A: (100, -0.5)}, **{cc.POSITIVE_CONTROL: (50, 0.1)})
    assert run.select(broken, _tests(**{cc.POSITIVE_CONTROL: 0.9}), 30)[0] == "readout_not_shown_to_work"

    supported = _table(**{**floor, run.ARM_A: (100, -0.5), run.ARM_C: (100, 0.0)}, **working)
    chosen = run.select(supported, _tests(**{cc.POSITIVE_CONTROL: 1e-6, run.ARM_A: 1e-4, run.ARM_C: 0.5}), 30)
    assert chosen[0] == "chain_supported"

    proximity = _table(**{**floor, run.ARM_A: (100, 0.01), run.ARM_C: (100, -0.4)}, **working)
    assert (
        run.select(proximity, _tests(**{cc.POSITIVE_CONTROL: 1e-6, run.ARM_A: 0.6, run.ARM_C: 1e-3}), 30)[0]
        == "proximity_explains_it"
    )

    negative = _table(**{**floor, run.ARM_A: (100, 0.01), run.ARM_C: (100, 0.02)}, **working)
    assert (
        run.select(negative, _tests(**{cc.POSITIVE_CONTROL: 1e-6, run.ARM_A: 0.6, run.ARM_C: 0.7}), 30)[0]
        == "no_difference_detected"
    )


def test_a_p_of_exactly_the_alpha_does_not_pass_and_a_median_tie_does_not_either():
    run = _script("closure_cell_run")
    table = _table(**{run.ARM_D: (200, 0.0), run.ARM_A: (100, -0.5), cc.POSITIVE_CONTROL: (50, -1.0)})
    at_alpha = run.select(table, _tests(**{cc.POSITIVE_CONTROL: 1e-6, run.ARM_A: 0.05}), 30)
    assert at_alpha[0] == "no_difference_detected"
    tied = _table(**{run.ARM_D: (200, 0.0), run.ARM_A: (100, 0.0), cc.POSITIVE_CONTROL: (50, -1.0)})
    assert (
        run.select(tied, _tests(**{cc.POSITIVE_CONTROL: 1e-6, run.ARM_A: 1e-8}), 30)[0]
        == "no_difference_detected"
    )


def test_the_reading_is_the_registered_wording_with_only_its_numbers_filled_in():
    run = _script("closure_cell_run")
    register = _script("closure_cell_register")
    table = _table(
        **{
            run.ARM_D: (200, 0.011),
            run.ARM_A: (100, -0.482),
            run.ARM_C: (90, -0.02),
            cc.POSITIVE_CONTROL: (50, -1.2),
        }
    )
    tests = _tests(**{cc.POSITIVE_CONTROL: 1e-9, run.ARM_A: 0.0004, run.ARM_C: 0.4})
    wording = register.OUTCOMES["chain_supported"]
    reading = run.fill(wording, table, tests, register.MIN_PER_ARM)
    assert "{" not in reading
    assert "100" in reading
    assert "200" in reading
    assert "-0.482" in reading
    # every word of the template that is not a placeholder survives verbatim
    for fragment in wording.replace("{", "|{").replace("}", "}|").split("|"):
        if fragment and not fragment.startswith("{"):
            assert fragment in reading


# ---- the local store, where it is present --------------------------------------------------------


@pytest.mark.needs_local_data("data/knowledge/closure_cell/ENCFF630CYG.bed.gz", how=GUIDE_BED_HOW)
def test_the_library_design_on_disk_is_the_file_encode_published_and_carries_no_count():
    assert hashlib.md5(cc.GUIDE_BED.read_bytes()).hexdigest() == cc.GUIDE_BED_MD5  # noqa: S324
    guides = cc.guides_by_chrom()
    assert guides, "the design BED is on disk but parsed to nothing"
    assert "chr1" in guides
    first = guides["chr1"][0]
    assert first[1] - first[0] == 3  # a 3 bp target interval, as the quantification files carry
    # six columns: position, name, score and strand. No count field exists in it to read.
    with gzip.open(cc.GUIDE_BED, "rt") as fh:
        assert len(fh.readline().rstrip("\n").split("\t")) == 6
