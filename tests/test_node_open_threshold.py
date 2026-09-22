# SPDX-License-Identifier: AGPL-3.0-or-later
"""The open-node call has one definition, and these pin what it is and what it is not.

`genomeos/attribution/candidates.py` re-derived `reader.py`'s median split line for line until
2026-09-22. The duplication was harmless - nothing in the repository reads the two fields it feeds
- but one rule with two definitions can be changed in one place and not the other, so the rule now
lives in `genomeos.genome.reader.node_open_threshold` and both sites call it. See
docs/NODES-READER-WRITER.md, "The second copy of that threshold, in the attribution candidates".
"""

from __future__ import annotations

import pytest

from genomeos.attribution import candidates
from genomeos.genome import reader


def _literal(densities: list[float]) -> float:
    """The expression that stood at reader.py:197-198 and candidates.py:401-404 before the merge."""
    d = sorted(densities)
    return max(1.0, d[len(d) // 2] if d else 0.0)


@pytest.mark.parametrize(
    "densities",
    [
        [],
        [0.0],
        [0.0, 0.0, 0.0],
        [0.5, 0.7, 0.9],  # every node under the floor: the floor, not the median, decides
        [1.0, 2.0, 3.0, 4.0],
        [float(i) for i in range(1, 1002)],
        [i * 0.01 for i in range(1, 500)],
    ],
)
def test_shared_threshold_reproduces_the_expression_it_replaced(densities):
    assert reader.node_open_threshold(densities) == _literal(densities)


def test_neither_module_grows_its_own_copy_back():
    """The names are imported, not re-derived, and the literal that was duplicated is gone."""
    assert candidates.node_open_threshold is reader.node_open_threshold
    for path in (candidates.__file__, reader.__file__):
        with open(path) as fh:
            body = fh.read()
        assert "max(1.0, median" not in body, f"the literal median split is back in {path}"


def test_the_call_is_a_within_biosample_rank_and_is_depth_invariant():
    """The property that makes this copy legitimate where the reader's `nodes_open` count was not:
    membership above the median survives scaling every density in a biosample, so a shallow
    experiment and a deep one rank the same node the same way."""
    densities = [float(i) for i in range(1, 1001)]
    base = reader.node_open_threshold(densities)
    open_ids = {i for i, d in enumerate(densities) if d >= base}
    for depth in (0.1, 2.0, 6.81):  # 6.81 is the observed DNase depth span over the thirteen
        scaled = [d * depth for d in densities]
        thr = reader.node_open_threshold(scaled)
        assert {i for i, d in enumerate(scaled) if d >= thr} == open_ids


def test_the_call_is_half_the_nodes_and_so_can_never_say_a_node_is_shut():
    """The other half of the same property, and the reason `nodes_open` was withdrawn as a count:
    a median split returns half whatever the depth, so the closed side carries no information."""
    for scale in (1, 10, 100):
        densities = [i * scale for i in range(1, 1001)]
        thr = reader.node_open_threshold(densities)
        assert sum(1 for d in densities if d >= thr) == 500


def test_below_the_floor_the_call_stops_being_a_rank():
    """The one case where this copy would produce a between-biosample artefact, kept as a named
    case rather than tidied away: on disk the floor binds only on chrY, where the open share runs
    0.00 to 0.49 across the thirteen biosamples instead of 0.50."""
    shallow = [0.0] * 40 + [0.1, 0.4, 0.9]  # median 0.0, every node under the floor
    deep = [d * 20 for d in shallow]  # the same biosample sequenced deeper
    assert reader.node_open_threshold(shallow) == reader.NODE_OPEN_MIN_DENSITY
    n_shallow = sum(1 for d in shallow if d >= reader.node_open_threshold(shallow))
    n_deep = sum(1 for d in deep if d >= reader.node_open_threshold(deep))
    assert n_shallow == 0 and n_deep == 3  # depth alone moved the count: not a rank here


def test_the_basis_travels_with_the_record():
    """A consumer opening a syntax_candidates_*.json must learn from the file that node_open is a
    rank, because nothing in the repository reads the field and no doc would be consulted."""
    assert reader.NODE_OPEN_BASIS is candidates.NODE_OPEN_BASIS
    for phrase in ("rank within one cell", "not an absolute call", "median", "never that the node is shut"):
        assert phrase in reader.NODE_OPEN_BASIS
