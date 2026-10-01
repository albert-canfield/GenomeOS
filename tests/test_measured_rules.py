"""The named overlap predicates: a discovery primitive, not an adopted attribution policy.

`measured.RECIPROCAL_HALF` is the production rule and decides evidence. `MIN_SIDE_HALF` is audit A's policy
(dd8c49c), offered as a predicate only: an overlap under it is a candidate and does not establish which
element produced a measured response. It is deliberately not reachable through `CrispriIndex.of()`, the
evidence interface, so one observation cannot become decisive evidence for several elements.

The implication these tests check is overlap inclusion. It is not preservation of unique attribution or of
verdict eligibility: a broader rule can keep every observation while resolving each to several elements.
"""

import random

from genomeos.attribution import ablation
from genomeos.attribution import measured as ms


def test_the_default_rule_is_the_production_one():
    assert ms.ATTACHMENT_RULES[ms.RECIPROCAL_HALF] is ms.measures
    assert ms.attaches(ms.RECIPROCAL_HALF, 1000, 1300, 1050, 1350) == ms.measures(1000, 1300, 1050, 1350)


def test_the_policy_attaches_a_wide_interval_the_production_rule_cannot():
    """The registered reason: a 300 bp registry element and a 1 kb tested interval differ by more than a
    factor of two, so the reciprocal rule can never pair them however they lie."""
    element, tested = (1000, 1300), (1000, 2000)
    assert ms.attaches(ms.RECIPROCAL_HALF, *element, *tested) is False
    assert ms.attaches(ms.MIN_SIDE_HALF, *element, *tested) is True


def test_the_production_rule_implies_the_policy_on_random_intervals():
    """Overlap inclusion only: half of both widths is at least half of the smaller, so the broader rule
    accepts every interval pair the production rule accepts. This says nothing about which element produced
    a response, and nothing about verdict eligibility."""
    rng = random.Random(20261001)
    checked = 0
    for _ in range(20000):
        a = rng.randrange(0, 2000)
        b = rng.randrange(0, 2000)
        element = (a, a + rng.randrange(1, 900))
        tested = (b, b + rng.randrange(1, 3000))
        if ms.attaches(ms.RECIPROCAL_HALF, *element, *tested):
            checked += 1
            assert ms.attaches(ms.MIN_SIDE_HALF, *element, *tested)
    assert checked > 100, "the property was not exercised"


def test_the_policy_uses_the_measured_layers_own_threshold_and_adds_no_second_one():
    assert ms.measures_min_side(0, 100, 50, 150) is (ms.min_side_overlap(0, 100, 50, 150) >= 0.5)
    assert ms.measures_min_side(0, 100, 50, 150, fraction=ms.RECIPROCAL_OVERLAP) is True


def test_min_side_overlap_is_zero_for_disjoint_or_empty_intervals():
    assert ms.min_side_overlap(0, 100, 200, 300) == 0.0
    assert ms.min_side_overlap(0, 0, 0, 100) == 0.0
    assert ms.min_side_overlap(0, 100, 50, 50) == 0.0


def test_min_side_overlap_is_one_when_the_smaller_interval_is_swallowed():
    assert ms.min_side_overlap(100, 200, 0, 1000) == 1.0


def test_an_unknown_rule_is_refused_rather_than_treated_as_the_default():
    """Silently falling back would let a result say it used a rule it did not."""
    try:
        ms.attaches("half_of_something", 0, 100, 0, 100)
    except ValueError as e:
        assert "unknown attachment rule" in str(e)
    else:
        raise AssertionError("an unknown rule must be refused")


def _pair(chrom: str, start: int, end: int) -> ms.CrispriPair:
    return ms.CrispriPair(
        chrom=chrom,
        start=start,
        end=end,
        gene="MYC",
        cell="K562",
        dataset="test",
        reference="test",
        regulated=True,
        significant=True,
        effect_size=-0.5,
        p_adjusted=0.01,
    )


def test_the_evidence_interface_does_not_expose_the_broader_rule():
    """`of()` decides evidence, so it takes no rule. A broader overlap would return one observation for
    several elements, making it decisive evidence for each; that needs a discovery interface with
    observation ids and whole candidate sets, which is not built."""
    import inspect

    assert "rule" not in inspect.signature(ablation.CrispriIndex.of).parameters


def test_the_evidence_interface_keeps_the_production_behaviour():
    """A wide tested interval the production rule cannot pair is not found, and no argument makes it so."""
    index = ablation.CrispriIndex([_pair("chr8", 1000, 2000)])
    assert index.of("chr8", 1000, 1300) == []
    assert ms.attaches(ms.MIN_SIDE_HALF, 1000, 1300, 1000, 2000) is True


def test_a_bad_rule_name_is_refused_before_any_interval_is_examined():
    """The empty cases: a name must not be answered with a plain "no overlap". Disjoint intervals, zero
    width and identical intervals all go through the same check."""
    for a_start, a_end, b_start, b_end in ((0, 10, 900, 1000), (5, 5, 0, 100), (0, 100, 0, 100)):
        try:
            ms.attaches("not_a_rule", a_start, a_end, b_start, b_end)
        except ValueError as e:
            assert "unknown attachment rule" in str(e)
        else:
            raise AssertionError("an unknown rule must be refused whatever the intervals")
