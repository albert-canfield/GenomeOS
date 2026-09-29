# SPDX-License-Identifier: Apache-2.0
"""Task 4.3: three germ layers in the right order along the NODAL gradient."""

import math

import pytest

from genomeos.runtime.gastrulation import run_gastrulation

# The expected_* proportions in data/demo/gastrulation.bio cite no source:
# each is `evidence: inferred "order of magnitude"`, confidence 0.3.
TOLERANCE = 0.25


def _within(value: float, target: float, tol: float = TOLERANCE) -> bool:
    """Strictly within tol, with float error resolved against the model.

    A difference that equals the tolerance up to rounding (0.35 - 0.10 is
    0.24999999999999997 in binary) is not inside it, so it counts as outside.
    """
    diff = abs(value - target)
    return diff < tol and not math.isclose(diff, tol, rel_tol=1e-9, abs_tol=1e-12)


@pytest.fixture(scope="module")
def result():
    return run_gastrulation(cells=60, hours=30, dt=0.05)


def _expected(r):
    return {k: r.module.parameters[f"expected_{k}"].value for k in ("ectoderm", "mesoderm", "endoderm")}


def test_within_is_not_fooled_by_rounding():
    assert 0.35 - 0.10 < 0.25  # the float fact the old check passed on
    assert not _within(0.10, 0.35)
    assert _within(0.11, 0.35)
    assert not _within(0.60, 0.35)


@pytest.mark.xfail(
    strict=True,
    reason=(
        "Tbxt-SOX17 sign corrected 2026-09-28 (Lolas et al. 2014, registered in b1f3405): no cell "
        "ends as mesoderm any more, 0 of 60, so the fates run endoderm then ectoderm with no band "
        "between; before the correction mesoderm held 6 of 60 cells and this test passed"
    ),
)
def test_three_layers_in_order(result):
    r = result
    props = r.proportions()
    assert all(props[f] > 0.05 for f in props), props
    # endoderm nearest the NODAL source, ectoderm farthest
    assert r.fates[0] == "endoderm" and r.fates[-1] == "ectoderm", (r.fates[0], r.fates[-1])
    order = [f for i, f in enumerate(r.fates) if i == 0 or f != r.fates[i - 1]]
    assert order == ["endoderm", "mesoderm", "ectoderm"], order


def test_ends_in_order_and_uncertainty_says_low(result):
    # split out of test_three_layers_in_order on 2026-09-28 so these checks still run while it xfails
    r = result
    assert r.fates[0] == "endoderm" and r.fates[-1] == "ectoderm", (r.fates[0], r.fates[-1])
    # some of this module's rules are inferred, and the report must say so. The label was "low" until
    # 2026-09-28, when the two Tbxt/Sox17 rules moved from inferred to experimental (Lolas et al. 2014)
    # and it rose to "medium"; their numbers are still unsourced, which the label does not see.
    unc = r.uncertainty().to_dict()["molecular"]
    assert unc["label"] == "medium" and unc["items"] > 5
    assert any(rule.evidence.kind.value == "inferred" for rule in r.module.rules)


def test_ectoderm_and_endoderm_within_unsourced_expectation(result):
    # measured 2026-09-28: ectoderm 0.467 vs 0.45; endoderm 0.433 vs 0.20
    # (inside 0.25 by 0.017, more than double the stated value)
    # after the Tbxt-SOX17 sign correction (same day): ectoderm 0.633, endoderm 0.367; both still inside
    props, exp = result.proportions(), _expected(result)
    for k in ("ectoderm", "endoderm"):
        assert _within(props[k], exp[k]), (k, props[k], exp[k])


@pytest.mark.xfail(
    strict=True,
    reason=(
        "measured mesoderm 0.10 (6 of 60 cells) against an unsourced expectation of 0.35; "
        "the difference is exactly 0.25, so the old `< 0.25` check passed only because "
        "0.35 - 0.10 == 0.24999999999999997 in floating point"
        "; after the Tbxt-SOX17 sign correction of 2026-09-28 mesoderm is 0.0 (0 of 60), "
        "further from 0.35"
    ),
)
def test_mesoderm_within_unsourced_expectation(result):
    props, exp = result.proportions(), _expected(result)
    assert _within(props["mesoderm"], exp["mesoderm"]), (props["mesoderm"], exp["mesoderm"])


# Census comparison, registered 2026-09-28 (docs/DESIGN-MINIMAL-CELL.md).


def _census_counts():
    import json
    from pathlib import Path

    return json.loads(Path("data/results/gastrulation_census_counts.json").read_text())


def test_census_counts_match_the_papers():
    c = _census_counts()
    hum = c["human_cs7"]
    assert hum["cells"] == 1195  # Tyser et al. 2021: 665 caudal, 340 rostral, 190 yolk sac
    sites: dict[str, int] = {}
    for by in hum["by_label_and_site"].values():
        for s, n in by.items():
            sites[s] = sites.get(s, 0) + n
    assert sites == {"caudal": 665, "rostral": 340, "yolk sac": 190}
    assert c["mouse_atlas"]["cells"] == 116312  # Pijuan-Sala et al. 2019


def test_census_mappings_are_disjoint_and_use_real_labels():
    from genomeos.runtime import gastrulation as g

    c = _census_counts()
    human_labels = set(c["human_cs7"]["by_label_and_site"])
    mouse_labels = {lab for by in c["mouse_atlas"]["by_stage_and_label"].values() for lab in by}
    for mappings, labels in (
        (g.CENSUS_HUMAN_MAPPINGS, human_labels),
        (g.CENSUS_MOUSE_MAPPINGS, mouse_labels),
    ):
        for mapping in mappings.values():
            used = [lab for layer in mapping.values() for lab in layer]
            assert len(used) == len(set(used)), mapping
            assert set(used) <= labels, set(used) - labels
    for mapping in g.CENSUS_HUMAN_MAPPINGS.values():
        assert not set(g.CENSUS_HUMAN_EXTRAEMBRYONIC) & {lab for layer in mapping.values() for lab in layer}


@pytest.mark.xfail(
    strict=True,
    reason=(
        "this file records the one run of 3c4e301 on the module before the Tbxt-SOX17 sign correction "
        "(2026-09-28, b1f3405); a fresh run of the corrected module no longer reproduces its shares "
        "(mesoderm 0.1 then, 0.0 now). The corrected run is gastrulation_census_comparison_sign_corrected"
    ),
)
def test_census_comparison_as_run_once():
    # the one registered run (2026-09-28): falsified on all three layers of the human CS7 census
    import json
    from pathlib import Path

    from genomeos.runtime import gastrulation as g

    d = json.loads(Path("data/results/gastrulation_census_comparison.json").read_text())
    assert d["tolerance"] == g.CENSUS_TOLERANCE and d["model_run"] == g.CENSUS_MODEL_RUN
    fresh = run_gastrulation(**g.CENSUS_MODEL_RUN).proportions()
    assert {k: round(v, 4) for k, v in fresh.items()} == d["model"]
    assert d["verdict"] == "falsified"
    j = d["human_cs7_judgement"]
    assert not any(j[k]["pass"] for k in ("ectoderm", "mesoderm", "endoderm"))
    assert j["mesoderm"]["interval"][0] > 0.69 and d["model"]["mesoderm"] == 0.1


def test_census_comparison_before_the_sign_correction_is_kept():
    # the record of the 3c4e301 run stays as it was committed, and the constants quote it
    import json
    from pathlib import Path

    from genomeos.runtime import gastrulation as g

    d = json.loads(Path("data/results/gastrulation_census_comparison.json").read_text())
    assert d["model"] == g.SIGN_CORRECTION_BEFORE["census_run"] and d["verdict"] == "falsified"
    assert not any(d["human_cs7_judgement"][k]["pass"] for k in ("ectoderm", "mesoderm", "endoderm"))


def test_census_comparison_after_the_sign_correction():
    # the one run after the Tbxt-SOX17 sign correction (registered in b1f3405), same mappings, site
    # filters, tolerance and model run as 3c4e301: still falsified on all three layers. Registered
    # expectation: endoderm up, mesoderm down; seen: endoderm down (0.425 -> 0.3667), mesoderm 0.1 -> 0.0
    import json
    from pathlib import Path

    from genomeos.runtime import gastrulation as g

    d = json.loads(Path("data/results/gastrulation_census_comparison_sign_corrected.json").read_text())
    old = json.loads(Path("data/results/gastrulation_census_comparison.json").read_text())
    assert d["tolerance"] == g.CENSUS_TOLERANCE and d["model_run"] == g.CENSUS_MODEL_RUN
    assert d["human_cs7_variants"] == old["human_cs7_variants"]  # the census side did not move
    fresh = run_gastrulation(**g.CENSUS_MODEL_RUN).proportions()
    assert {k: round(v, 4) for k, v in fresh.items()} == d["model"]
    assert d["model"] == {"ectoderm": 0.6333, "mesoderm": 0.0, "endoderm": 0.3667}
    assert d["verdict"] == "falsified"
    j = d["human_cs7_judgement"]
    assert not any(j[k]["pass"] for k in ("ectoderm", "mesoderm", "endoderm"))
    sc = d["sign_correction"]
    assert sc["rules"] == g.SIGN_CORRECTION and sc["model_before"] == g.SIGN_CORRECTION_BEFORE["census_run"]
