# SPDX-License-Identifier: AGPL-3.0-or-later
"""The diagnosis lane's own logic: the four classes, the frozen ranking, and the committed counts."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

from genomeos.attribution import measured as ms

ROOT = Path(__file__).resolve().parents[1]


def _module(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pc = _module("pc198", "scripts/placement_cause_198.py")
pa = pc.pa


def _row(locus: int, cause: str, start: int = 0, end: int = 500):
    pair = ms.CrispriPair(
        chrom="chr1",
        start=start,
        end=end,
        gene="G",
        cell="K562",
        dataset="D",
        reference="R",
        regulated=True,
        significant=True,
        effect_size=-0.3,
        p_adjusted=0.01,
    )
    return (locus, locus, pair, cause)


# --- the classification is imported, not invented --------------------------------------------------
def test_the_causes_and_their_order_are_audit_a_s_own() -> None:
    """The lane imports audit A rather than restating it: the module it reads them from is that file."""
    assert Path(pa.__file__).resolve() == (ROOT / "scripts/placement_audit.py").resolve()
    assert pa.CAUSES[0] == pa.PLACED  # the order is the cascade's, first match wins
    assert set(pa.CAUSE_MEANING) == set(pa.CAUSES)
    assert pa.EXCLUSION[0] == "compiled_under_other_coordinates"


def test_every_cause_but_placed_has_exactly_one_class() -> None:
    assert set(pc.CLASS_OF) == set(pa.CAUSES) - {pa.PLACED}
    assert set(pc.CLASS_OF.values()) <= set(pc.CLASSES)
    assert len(pc.CLASSES) == len(set(pc.CLASSES)) == 4
    assert set(pc.CLASS_MEANING) == set(pc.CLASSES)


def test_class_four_is_named_and_is_not_a_catch_all() -> None:
    """Only the two causes that are neither an absence nor a width failure may land in class 4, and
    the class says which they are."""
    fourth = {c for c, k in pc.CLASS_OF.items() if k == pc.CLASSES[3]}
    assert fourth == {pa.PARTIAL, pa.ZERO_WIDTH}
    assert pa.PARTIAL in pc.CLASS_MEANING[pc.CLASSES[3]]
    assert pa.ZERO_WIDTH in pc.CLASS_MEANING[pc.CLASSES[3]]


def test_width_failures_and_absence_are_not_the_same_class() -> None:
    assert pc.CLASS_OF[pa.UNREACHABLE] == pc.CLASS_OF[pa.WIDTH_MISMATCH] == pc.CLASSES[1]
    assert pc.CLASS_OF[pa.ABSENT] == pc.CLASSES[0]
    assert pc.CLASS_OF[pa.REGISTRY_WOULD] == pc.CLASSES[2]


# --- one cause per locus, by the frozen order ------------------------------------------------------
def test_a_locus_takes_the_earliest_cause_of_the_frozen_order_any_pair_reaches() -> None:
    order = {c: i for i, c in enumerate(pa.CAUSES)}
    rows = [_row(1, pa.ABSENT), _row(1, pa.REGISTRY_WOULD), _row(1, pa.UNREACHABLE), _row(2, pa.PARTIAL)]
    assert pc.by_cause(rows, order) == {1: pa.REGISTRY_WOULD, 2: pa.PARTIAL}


def test_the_sensitivity_arm_takes_the_last_cause_instead() -> None:
    order = {c: i for i, c in enumerate(pa.CAUSES)}
    rows = [_row(1, pa.ABSENT), _row(1, pa.REGISTRY_WOULD), _row(1, pa.UNREACHABLE)]
    assert pc.by_cause(rows, order, last=True) == {1: pa.ABSENT}


def test_the_classes_sum_to_the_number_of_loci() -> None:
    causes = {1: pa.ABSENT, 2: pa.UNREACHABLE, 3: pa.WIDTH_MISMATCH, 4: pa.REGISTRY_WOULD, 5: pa.PARTIAL}
    got = pc.classes_of(causes)
    assert sum(got.values()) == len(causes)
    assert got == {
        pc.CLASSES[0]: 1,
        pc.CLASSES[1]: 2,
        pc.CLASSES[2]: 1,
        pc.CLASSES[3]: 1,
    }


def test_widths_report_the_bar_the_reciprocal_rule_cannot_pass() -> None:
    rows = [_row(1, pa.ABSENT, 0, 100), _row(1, pa.ABSENT, 0, 800), _row(1, pa.ABSENT, 0, 50)]
    got = pc.widths([p for _i, _g, p, _c in rows])
    assert (got["min"], got["max"]) == (50, 800)
    assert got["wider_than_700"] == 1
    assert got["narrower_than_75"] == 1


# --- the named alternative is the imported one, and nothing is adopted ------------------------------
def test_the_named_alternative_is_the_measured_layer_s_own_predicate() -> None:
    assert ms.ATTACHMENT_RULES[ms.MIN_SIDE_HALF] is ms.measures_min_side
    # the production rule is untouched by this lane
    assert ms.ATTACHMENT_RULES[ms.RECIPROCAL_HALF] is ms.measures
    assert ms.RECIPROCAL_OVERLAP == 0.5


# --- the committed counts --------------------------------------------------------------------------
def _committed() -> dict:
    return json.loads((ROOT / "data/results/placement_cause_198.json").read_text())


def test_the_committed_classes_sum_to_198_and_every_locus_has_one() -> None:
    d = _committed()
    classes = d["answer"]["classes"]
    assert sum(classes.values()) == d["answer"]["classes_sum"] == 198
    assert sum(d["answer"]["locus_causes"].values()) == 198
    assert d["population"]["independent_loci"] == 198
    assert d["population"]["unclassifiable_from_what_is_on_disk"] == 0


def test_the_committed_locus_causes_map_onto_the_committed_classes() -> None:
    d = _committed()
    folded: dict[str, int] = {}
    for cause, n in d["answer"]["locus_causes"].items():
        folded[pc.CLASS_OF[cause]] = folded.get(pc.CLASS_OF[cause], 0) + n
    assert folded == {k: v for k, v in d["answer"]["classes"].items() if v}


def test_the_committed_gates_reproduced_both_earlier_results() -> None:
    d = _committed()
    assert d["gates"]["equal_to_committed"] is True
    assert d["gates"]["equal_to_committed_audit_a"] is True
    assert d["gates"]["compiled_element_sets"]["pairs_whose_cause_differs_between_them"] == 0
    assert d["gates"]["response_map_coverage_reproduced"] == {
        "measured_perturbation_pairs": 14734,
        "independent_loci": 473,
        "no_compiled_link_pairs": 9909,
        "no_compiled_link_loci": 198,
    }


def test_the_committed_result_changed_no_rule_and_spent_nothing() -> None:
    d = _committed()
    assert d["alphagenome_requests"] == 0
    assert d["downloads"] == 0
    assert d["per_element_response_cache_opened"] is False
    assert d["cost"]["requests"] == 0
    assert "not_adopted" in d["named_alternative_descriptively"]


def test_the_committed_manifest_is_complete_and_records_each_input_on_its_own_path() -> None:
    m = _committed()["result_manifest"]
    assert m["complete"] is True
    paths = [i["path"] for i in m["inputs"]]
    assert len(paths) == len(set(paths))
    assert not any("(" in p for p in paths)  # a group entry, never used by this lane
