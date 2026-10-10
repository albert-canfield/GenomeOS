# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-strengthform: the strength-dependent form, its algebra, and the engine it must not change.

THE ENGINE IS NOT THIS LANE'S TO EDIT. `genomeos/runtime/grn.py` is an Apache-2.0 BioLang engine
file, and the form under test is already in it as the opt-in rule `competition`. The first tests
here pin its licence header, its default rule, its two rule tuples and every existing form's
parameter-free identity, so a change to any of them fails this lane's own acceptance.

THE REST IS THE LANE'S ALGEBRA, planted rather than asserted: a single must come back exactly, two
equal maximal singles must land strictly between the stronger single and the additive expectation,
each fixed rule must give the closed form its own algebra gives, an interval rule must cost nothing
inside its interval, and the reading must render each of its three words from the intervals alone.
"""

from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path

import pytest

from genomeos.attribution import combine2 as c2
from genomeos.attribution import strengthform as sf
from genomeos.runtime import grn

ENGINE = Path("genomeos/runtime/grn.py")
LOCAL_HOW = "uv run python scripts/strengthform_fetch.py"


# --- the engine, unchanged -------------------------------------------------------------------


def test_the_engine_file_keeps_its_apache_header():
    head = ENGINE.read_text().splitlines()[:2]
    assert head[0] == "# SPDX-License-Identifier: Apache-2.0"
    assert "BioLang engine" in head[1]


def test_the_engine_default_rule_and_both_rule_tuples_are_unchanged():
    assert grn.ACTIVATOR_COMBINATION == "mean"
    assert grn.ACTIVATOR_RULES == ("mean", "sum_capped", "max", "or")
    assert grn.OPT_IN_ACTIVATOR_RULES == ("competition", "multiplicative")
    assert sf.FORM in grn.OPT_IN_ACTIVATOR_RULES


def test_every_existing_form_keeps_its_identity():
    terms = [0.2, 0.5]
    assert grn.combine_activators(terms, "mean") == pytest.approx(0.35)
    assert grn.combine_activators(terms, "sum_capped") == pytest.approx(0.7)
    assert grn.combine_activators(terms, "max") == pytest.approx(0.5)
    assert grn.combine_activators(terms, "or") == pytest.approx(1 - 0.8 * 0.5)
    assert grn.combine_activators(terms, "competition", 1.0) == pytest.approx(0.7 / 1.7)
    assert grn.combine_activators(terms, "multiplicative", gain=1.0) == pytest.approx(math.exp(0.7))
    assert grn.combine_activators([], "competition", 1.0) == 1.0  # no activator: unchanged


def test_the_lane_declares_the_same_theta_as_7d():
    assert sf.THETA == c2.THETA == 0.5
    assert sf.HALF_SATURATION == 1.0


# --- the form's algebra ----------------------------------------------------------------------


@pytest.mark.parametrize("bmax", [0.25, 1.0, 3.0])
@pytest.mark.parametrize("share", [0.1, 0.5, 1.0])
def test_a_single_comes_back_exactly(bmax, share):
    b = bmax * share
    assert sf.predicted_boost(sf.FORM, b, 1e-12, bmax) == pytest.approx(b, abs=1e-9)


def test_the_strongest_single_has_term_one_at_the_declared_theta():
    t1, t2 = sf.terms(2.0, 2.0, 2.0)
    assert t1 == pytest.approx(1.0)
    assert t2 == pytest.approx(1.0)
    assert grn.combine_activators([t1], sf.FORM, sf.HALF_SATURATION) == pytest.approx(sf.THETA)


def test_two_equal_maximal_singles_lie_between_the_stronger_single_and_additive():
    b = 1.0
    own = sf.predicted_boost(sf.FORM, b, b, b)
    lower, upper = sf.predicted_boost("additive_envelope", b, b, b)
    assert lower == pytest.approx(b)
    assert upper == pytest.approx(math.log2(2 * 2**b - 1))
    assert lower < own < upper


def test_the_form_is_strength_dependent_and_not_a_fixed_function_of_the_terms():
    """Two pairs with the same singles but different promoter maxima predict differently."""
    weak = sf.predicted_boost(sf.FORM, 0.2, 0.2, 3.0)
    strong = sf.predicted_boost(sf.FORM, 0.2, 0.2, 0.2)
    additive = math.log2(2 * 2**0.2 - 1)
    assert weak > strong
    assert weak == pytest.approx(additive, abs=0.02)  # weak engagement adds
    assert strong < additive - 0.02  # strong engagement saturates


def test_each_fixed_rule_gives_its_own_closed_form():
    b1, b2 = 0.5, 1.5
    assert sf.predicted_boost("no_rule", b1, b2, 2.0) == 0.0
    assert sf.predicted_boost("mean", b1, b2, 2.0) == pytest.approx(math.log2((2**b1 + 2**b2) / 2))
    assert sf.predicted_boost("max", b1, b2, 2.0) == pytest.approx(b2)
    assert sf.predicted_boost("multiplicative", b1, b2, 2.0) == pytest.approx(b1 + b2)
    lower, upper = sf.predicted_boost("additive_envelope", b1, b2, 2.0)
    assert lower == pytest.approx(b2)
    assert upper == pytest.approx(math.log2(2**b1 + 2**b2 - 1))


def test_an_unknown_form_and_an_impossible_occupancy_are_refused():
    with pytest.raises(ValueError):
        sf.predicted_boost("sum_capped", 0.5, 0.5, 1.0)
    with pytest.raises(ValueError):
        sf.terms(0.5, 0.5, 0.0)  # a promoter with no activating single
    with pytest.raises(ValueError):
        sf.terms(2.0, 0.5, 1.0)  # a single above the promoter's maximum


def test_an_interval_rule_costs_nothing_inside_its_interval():
    assert sf.squared_error(1.2, (1.0, 1.5)) == 0.0
    assert sf.squared_error(0.8, (1.0, 1.5)) == pytest.approx(0.04)
    assert sf.squared_error(2.0, (1.0, 1.5)) == pytest.approx(0.25)
    assert sf.squared_error(2.0, 1.0) == pytest.approx(1.0)


# --- the reading -----------------------------------------------------------------------------


def _intervals(delta: float, lo: float, hi: float) -> dict[str, dict[str, float]]:
    row = {"delta_mse": delta, "ci_low": lo, "ci_high": hi, "draws": 10}
    return {name: dict(row) for name in sf.COMPARATORS}


def test_each_reading_word_comes_from_the_interval_alone():
    assert sf.word({"ci_low": 0.1, "ci_high": 0.2}) == "better"
    assert sf.word({"ci_low": -0.2, "ci_high": -0.1}) == "worse"
    assert sf.word({"ci_low": -0.1, "ci_high": 0.2}) == "no difference detected"


def test_the_reading_renders_better_and_is_then_not_falsified():
    out = sf.reading(_intervals(0.3, 0.1, 0.5), own_mse=1.0, n=10, clusters=3, promoters=2)
    assert out["better_than_every_comparator"] is True
    assert out["falsified"] is False
    assert sf.VERDICT_BETTER in out["sentence"]
    assert sf.FALSIFIED_NO in out["sentence"]


def test_one_comparator_it_does_not_beat_falsifies_the_claim():
    intervals = _intervals(0.3, 0.1, 0.5)
    intervals["mean"] = {"delta_mse": -0.2, "ci_low": -0.3, "ci_high": -0.1, "draws": 10}
    out = sf.reading(intervals, own_mse=1.0, n=10, clusters=3, promoters=2)
    assert out["words"]["mean"] == "worse"
    assert out["falsified"] is True
    assert sf.VERDICT_NOT_BETTER in out["sentence"]
    assert sf.FALSIFIED_YES in out["sentence"]
    assert "mean" in out["sentence"]


def test_the_bootstrap_is_deterministic_and_resamples_whole_clusters():
    clusters = [("a", "b"), ("a", "b"), ("c", "d"), ("c", "d"), ("e", "f")]
    errors = {
        sf.FORM: [1.0, 1.0, 2.0, 2.0, 3.0],
        **{name: [2.0, 2.0, 3.0, 3.0, 4.0] for name in sf.COMPARATORS},
    }
    first = sf.cluster_bootstrap(clusters, errors, draws=50)
    again = sf.cluster_bootstrap(clusters, errors, draws=50)
    assert first == again
    for row in first.values():
        assert row["delta_mse"] == pytest.approx(1.0)
        assert row["ci_low"] == pytest.approx(1.0)
        assert row["ci_high"] == pytest.approx(1.0)


# --- the table reader ------------------------------------------------------------------------

COLUMNS = ("class", "Promoter", "frag1", "strand1", "frag2", "strand2", "Baseline", "Boost1", "Boost2")


def _planted(path: Path, rows: list[dict[str, str]]) -> Path:
    with gzip.open(path, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[*COLUMNS, "id", sf.OUTCOME_COLUMN], delimiter="\t")
        writer.writeheader()
        for i, row in enumerate(rows):
            writer.writerow({"id": f"row{i}", **row})
    return path


def _row(**over: str) -> dict[str, str]:
    base = {
        "class": "EE",
        "Promoter": "P1",
        "frag1": "E1",
        "strand1": "+",
        "frag2": "E2",
        "strand2": "+",
        "Baseline": "0.5",
        "Boost1": "1.0",
        "Boost2": "0.5",
        sf.OUTCOME_COLUMN: "1.0",
    }
    return {**base, **over}


def test_the_reader_keeps_only_informative_combinations_and_counts_the_rest(tmp_path):
    path = _planted(
        tmp_path / "t.tsv.gz",
        [
            _row(),
            _row(**{"class": "ER", "frag2": "RandomE1"}),
            _row(Boost1="NA"),
            _row(Boost2="-0.5", frag2="E3"),
            _row(**{sf.OUTCOME_COLUMN: "0"}),
            _row(Baseline="0"),
        ],
    )
    table = sf.read_table(path, read_outcome=True)
    counts = table["counts"]
    assert counts["rows"] == 6
    assert counts["class_not_EE"] == 1
    assert counts["single_missing"] == 1
    assert counts["single_not_above_baseline"] == 1
    assert counts["activity_not_positive"] == 1
    assert counts["no_baseline"] == 1
    assert counts["informative"] == 1
    kept = table["combinations"][0]
    assert kept.cluster == ("E1", "E2")
    assert kept.observed_boost == pytest.approx(1.0)
    assert table["promoter_maxima"]["P1"] == pytest.approx(1.0)


def test_the_reader_opens_no_outcome_when_it_is_not_asked_to(tmp_path):
    path = _planted(tmp_path / "t.tsv.gz", [_row(), _row(**{sf.OUTCOME_COLUMN: "0"})])
    table = sf.read_table(path, read_outcome=False)
    assert table["counts"]["informative"] == 2  # the non-positive activity was never looked at
    assert all(c.activity is None for c in table["combinations"])
    with pytest.raises(ValueError):
        _ = table["combinations"][0].observed_boost


# --- the real table, where it is on this machine ---------------------------------------------


@pytest.mark.needs_local_data(str(sf.STORE), how=LOCAL_HOW)
def test_the_pinned_table_is_the_one_on_disk():
    import hashlib

    h = hashlib.sha256()
    with sf.TABLE.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    assert h.hexdigest() == sf.TABLE_SHA256
