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


def test_three_layers_in_order(result):
    r = result
    props = r.proportions()
    assert all(props[f] > 0.05 for f in props), props
    # endoderm nearest the NODAL source, ectoderm farthest
    assert r.fates[0] == "endoderm" and r.fates[-1] == "ectoderm", (r.fates[0], r.fates[-1])
    order = [f for i, f in enumerate(r.fates) if i == 0 or f != r.fates[i - 1]]
    assert order == ["endoderm", "mesoderm", "ectoderm"], order
    # half of this module's rules are inferred, and the report must say so
    unc = r.uncertainty().to_dict()["molecular"]
    assert unc["label"] == "low" and unc["items"] > 5
    assert any(rule.evidence.kind.value == "inferred" for rule in r.module.rules)


def test_ectoderm_and_endoderm_within_unsourced_expectation(result):
    # measured 2026-09-28: ectoderm 0.467 vs 0.45; endoderm 0.433 vs 0.20
    # (inside 0.25 by 0.017, more than double the stated value)
    props, exp = result.proportions(), _expected(result)
    for k in ("ectoderm", "endoderm"):
        assert _within(props[k], exp[k]), (k, props[k], exp[k])


@pytest.mark.xfail(
    strict=True,
    reason=(
        "measured mesoderm 0.10 (6 of 60 cells) against an unsourced expectation of 0.35; "
        "the difference is exactly 0.25, so the old `< 0.25` check passed only because "
        "0.35 - 0.10 == 0.24999999999999997 in floating point"
    ),
)
def test_mesoderm_within_unsourced_expectation(result):
    props, exp = result.proportions(), _expected(result)
    assert _within(props["mesoderm"], exp["mesoderm"]), (props["mesoderm"], exp["mesoderm"])
