# SPDX-License-Identifier: AGPL-3.0-or-later
"""`!=` in a `when` clause is refused (lane-ne, 2026-09-28).

The parser used to turn `k != v` into the value `!=v`, which no matcher implements, so the clause
matched only the literal string "!=v" on every block. No program used it
(data/results/when_census.json: 0 of 456,157 clauses). Supporting it would mean deciding what a
missing key means, which no program needs, so the parser refuses it and names the alternatives.
"""

from __future__ import annotations

import pytest

from genomeos.lang import parse
from genomeos.lang.parser import BioLangError

GENES = (
    'gene A { symbol: A; evidence: curated "x"; confidence: 0.9 }\n'
    'gene B { symbol: B; evidence: curated "x"; confidence: 0.9 }\n'
)

BLOCKS = {
    "decision": "decision d {{ action: divide; when: {w}; daughters: X, Y; after: 1 min }}\n",
    "timer": "timer t {{ duration: 10 min; when: {w} }}\n",
    "rule": 'rule A activates B {{ when: {w}; evidence: curated "x"; confidence: 0.5 }}\n',
    "event": "event e {{ rate: 2 /yr; when: {w}; effect: damage -= 1 }}\n",
}

SPELLINGS = ["cell_type != K562", "cell_type = !=K562", "lineage = AB, cell_type != K562"]


@pytest.mark.parametrize("kind", sorted(BLOCKS))
@pytest.mark.parametrize("when", SPELLINGS)
def test_not_equal_is_refused_with_the_clause_and_the_alternatives(kind, when):
    src = "module t\n" + GENES + BLOCKS[kind].format(w=when)
    with pytest.raises(BioLangError) as err:
        parse(src)
    msg = str(err.value)
    assert "cell_type" in msg and "K562" in msg
    assert "`!=` is not supported" in msg
    for alternative in ("a|b", "absent", ">=n"):
        assert alternative in msg


@pytest.mark.parametrize("kind", sorted(BLOCKS))
def test_the_other_constructs_still_parse_on_every_block(kind):
    when = "cell_type = K562|HepG2, host = absent, generation >= 3, stage <= 2"
    parse("module t\n" + GENES + BLOCKS[kind].format(w=when))
