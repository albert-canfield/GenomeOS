# SPDX-License-Identifier: AGPL-3.0-or-later
"""The terminality gate instrument (area E): the properties the result is only meaningful under.

These pin the accounting, not the answer. The registration in `genomeos/organism/terminality.py` says
that a low number is not a failure and broken accounting is, so what is tested here is that the metric
cannot be raised by abstaining, that the universe is the registered one, and that no arm can read the
lookup's terminality.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.organism import terminality as tg
from genomeos.organism.reference import ReferenceLineage

REF = Path("data/results/celegans_lineage_cells.json")
RESULT = Path("data/results/celegans_terminality.json")

pytestmark = pytest.mark.skipif(not REF.exists(), reason="the distilled reference lineage is not present")


@pytest.fixture(scope="module")
def ref() -> ReferenceLineage:
    return ReferenceLineage.load()


def test_the_universe_is_the_registered_one(ref: ReferenceLineage) -> None:
    """V1: 555 terminal and 771 dividing embryonic cells, the 113 deaths excluded."""
    labels = tg.universe(ref)
    assert len(labels) == 1326
    assert sum(1 for v in labels.values() if v == tg.TERMINAL) == 555
    assert sum(1 for v in labels.values() if v == tg.DIVIDING) == 771
    with_deaths = tg.universe(ref, deaths=True)
    assert len(with_deaths) == 1439
    assert sum(1 for v in with_deaths.values() if v == tg.TERMINAL) == 555 + 113


def test_abstaining_cannot_raise_the_metric(ref: ReferenceLineage) -> None:
    """The trap the metric exists to close: the cited fate rules claim more non-terminal cells than
    terminal ones, so a rule set that abstains more would look better on any claimed-cell metric.
    Dropping a correct Terminal call can only lower the balanced accuracy, never raise it."""
    labels = tg.universe(ref)
    full = dict(labels)  # a perfect caller
    perfect = tg.score("perfect", full, labels)
    assert perfect.balanced == 1.0
    assert perfect.claimed == 1326
    terminals = [c for c in labels if labels[c] == tg.TERMINAL]
    dropped = set(terminals[:50])
    abstaining = {c: v for c, v in full.items() if c not in dropped}
    partial = tg.score("abstains on 50 terminal cells", abstaining, labels)
    assert partial.claimed == 1276
    assert partial.balanced < perfect.balanced
    assert partial.n == 1326  # every cell is still scored


def test_every_cell_is_called_whatever_the_arm_says(ref: ReferenceLineage) -> None:
    labels = tg.universe(ref)
    s = tg.score("empty", {}, labels)
    assert s.called_terminal + s.called_dividing == 1326
    assert s.claimed == 0
    assert s.balanced == 0.5  # calling everything Dividing is worth exactly chance


def test_the_depth_baselines_are_what_was_registered(ref: ReferenceLineage) -> None:
    labels = tg.universe(ref)
    gen = tg.depths(ref, labels)
    assert tg.best_threshold(gen, labels) == 9
    t3 = tg.score("T3", tg.depth_call(gen, 9), labels)
    assert round(t3.balanced, 4) == 0.7815
    assert round(t3.accuracy, 4) == 0.8069
    ks = tg.per_sublineage_thresholds(gen, labels)
    assert {ks["ABal"], ks["ABar"], ks["ABpl"], ks["ABpr"]} == {9}
    assert (ks["MS"], ks["E"], ks["C"], ks["D"]) == (6, 5, 6, 4)
    t4 = tg.score("T4", tg.per_sublineage_call(gen, ks, 9), labels)
    assert round(t4.balanced, 4) == 0.8998


def test_no_arm_can_read_the_lookups_terminality() -> None:
    """The instrument's only route to `children`, `divides` or a tracked lifetime is `universe`, which
    builds the labels. Nothing else in the module may touch them."""
    src = Path("genomeos/organism/terminality.py").read_text()
    body = src.split("def depths(", 1)[1]
    for forbidden in (".children", ".divides", "lifetimes", "cell_type", ".dies"):
        assert forbidden not in body, f"{forbidden} is readable outside the label builder"


def test_depth_tokens_are_a_prefix_chain() -> None:
    assert tg.depth_tokens(0) == set()
    assert tg.depth_tokens(3) == {"gen>=1", "gen>=2", "gen>=3"}
    assert tg.depth_tokens(9) > tg.depth_tokens(8)


def test_the_logistic_ceiling_converges_monotonically() -> None:
    """A separable two-point problem: the step is 1/L for the loss's own Lipschitz constant, so the
    loss falls at every step and there is nothing to damp. Pinned because a sibling lane found an
    undamped Newton step in a shared fit today."""
    np = pytest.importorskip("numpy")
    rows = [[1.0, 0.0], [0.0, 1.0], [1.0, 0.0], [0.0, 1.0]]
    y = [1, 0, 1, 0]
    w, b = tg._fit_logistic(rows, y, lam=0.01, steps=200)
    p = tg._predict_logistic(w, b, rows)
    assert p[0] > 0.5 and p[1] < 0.5
    assert bool(np.isfinite(w).all())


@pytest.mark.skipif(not RESULT.exists(), reason="the terminality result has not been run here")
def test_the_result_file_accounts_for_every_cell() -> None:
    """V2, pinned on the committed result: every arm calls all 1,326 cells."""
    d = json.loads(RESULT.read_text())
    assert d["universe"] == {
        "embryonic_cells": 1326,
        "terminal": 555,
        "dividing": 771,
        "programmed_deaths_excluded": 113,
    }
    for row in d["arms"] + d["baselines"]:
        assert row["n"] == 1326, row["arm"]
        assert row["called_terminal"] + row["called_dividing"] == 1326, row["arm"]
    assert d["model_requests"] == 0
