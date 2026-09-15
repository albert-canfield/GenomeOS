# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""`share`: a partition whose numbers are the published ones (BIOLANG-v0.3.md, populations).

A splitting `fraction` is a share of what is *left* when it runs, so sibling splits multiply and the
value of every split after the first depends on the order they are written in. That is how the human
body program came to carry numbers that its own evidence lines called shares of a germ layer while
being up to 140,000x away from one. `share` is a share of the population as it stood at the start of
the decision point — or of the instant, for recurring flows — so siblings are independent of each
other, and a program states the number the paper states.

The falsifier is the body program itself: restating its germ-layer splits as shares must not move the
solved steady state it already reproduces.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pytest

from genomeos.ir.model import to_minutes
from genomeos.lang import parse, parse_file
from genomeos.lang.parser import BioLangError
from genomeos.runtime.body import Body

ROOT = Path(__file__).resolve().parent.parent
BODY = ROOT / "data" / "organisms" / "human" / "body.bio"
YEARS20 = to_minutes(20, "yr")

POOL = """
module toy
import bio.std.development
cell_type A { parent: PostMitotic }
cell_type B { parent: PostMitotic }
cell_type C { parent: PostMitotic }
organism T { root: Pool; cell_type: Blastomere; resolution: populations; observe: count, fates }
stage S { from: 0 min; to: 10 min }
stage Split { from: 10 min }
decision to_A { action: differentiate; when: cell_type = Blastomere, stage = Split; to: A; HOW_A }
decision to_B { action: differentiate; when: cell_type = Blastomere, stage = Split; to: B; HOW_B }
decision to_C { action: differentiate; when: cell_type = Blastomere, stage = Split; to: C; HOW_C }
"""

FLOWS = """
module toy
import bio.std.development
cell_type A { parent: PostMitotic }
cell_type B { parent: PostMitotic }
organism T { root: Pool; cell_type: Blastomere; resolution: populations; observe: count }
stage S { from: 0 d }
decision out_A { action: differentiate; when: cell_type = Blastomere; to: A; HOW_A; after: 1 d }
decision out_B { action: differentiate; when: cell_type = Blastomere; to: B; HOW_B; after: 1 d }
"""


def census(source: str, until: float, **how: str) -> dict[str, float]:
    text = source
    for key, value in how.items():
        text = text.replace(key, value)
    body = Body(parse(text), seed=None).run(until=until)
    out: dict[str, float] = defaultdict(float)
    for c in body.alive_at(until):
        out[c.cell_type] += c.count
    return {k: v for k, v in out.items() if v > 1e-12}


# ---- what it means -------------------------------------------------------------------


def test_the_same_partition_written_both_ways_gives_the_same_split():
    """Sequential fractions and absolute shares describing one partition agree; only one of them
    states the numbers the partition is actually made of."""
    remainders = census(POOL, 20, HOW_A="fraction: 0.25", HOW_B="fraction: 0.3333", HOW_C="fraction: 1.0")
    shares = census(POOL, 20, HOW_A="share: 0.25", HOW_B="share: 0.25", HOW_C="share: 0.5")
    assert remainders.keys() == shares.keys()
    for kind, value in shares.items():
        assert remainders[kind] == pytest.approx(value, rel=1e-3)


def test_a_share_does_not_depend_on_where_it_is_written():
    """The order-dependence §7.3 removed for fates, one layer down: each decision takes what it says."""
    one = census(POOL, 20, HOW_A="share: 0.5", HOW_B="share: 0.25", HOW_C="share: 0.25")
    assert one == {"A": pytest.approx(0.5), "B": pytest.approx(0.25), "C": pytest.approx(0.25)}
    other = census(POOL, 20, HOW_A="share: 0.25", HOW_B="share: 0.25", HOW_C="share: 0.5")
    assert other == {"A": pytest.approx(0.25), "B": pytest.approx(0.25), "C": pytest.approx(0.5)}


def test_sibling_flows_at_one_instant_draw_from_one_base():
    """What haematopoiesis needs: outflows due together retain 1 - sum(s), not the product of
    (1 - f), so a generator no longer has to solve its growth term against the order the runtime
    happens to apply them in."""
    with_fractions = census(FLOWS, 1441, HOW_A="fraction: 0.1", HOW_B="fraction: 0.1")
    with_shares = census(FLOWS, 1441, HOW_A="share: 0.1", HOW_B="share: 0.1")
    assert with_fractions["Blastomere"] == pytest.approx(0.9 * 0.9)
    assert with_shares["Blastomere"] == pytest.approx(1.0 - 0.1 - 0.1)


def test_shares_that_ask_for_more_than_there_is_are_refused():
    """Absolute shares can over-subscribe, which is a program stating impossible facts. Refused with
    the numbers, never rescaled into something that runs."""
    with pytest.raises(ValueError, match="must not sum above 1"):
        census(POOL, 20, HOW_A="share: 0.6", HOW_B="share: 0.6", HOW_C="share: 0.1")


def test_the_grammar_refuses_what_it_cannot_mean():
    both = POOL.replace("HOW_A", "share: 0.5; fraction: 0.5").replace("HOW_B", "share: 0.2")
    with pytest.raises(BioLangError, match="not both"):
        parse(both.replace("HOW_C", "share: 0.3"))
    with pytest.raises(BioLangError, match="only a differentiate decision"):
        parse("module t\ndecision d { action: divide; when: cell = P0; share: 0.5 }\n")
    with pytest.raises(BioLangError, match="share must be within"):
        parse(
            POOL.replace("HOW_A", "share: 1.5").replace("HOW_B", "share: 0.2").replace("HOW_C", "share: 0.3")
        )


# ---- the falsifier: a solved steady state must not move ------------------------------


def _germ_layer_splits(module) -> int:
    """Restate each germ-layer split as the share of the layer its fraction was derived from."""
    groups: dict[str, list] = defaultdict(list)
    for d in module.decisions:
        if d.action == "differentiate" and "_to_" in d.id:
            groups[d.id.split("_to_")[0]].append(d)
    restated = 0
    for members in groups.values():
        left = 1.0
        for d in members:
            share = left * d.fraction
            left -= share
            d.share, d.fraction = share, 1.0
            restated += 1
    return restated


def test_restating_the_body_as_shares_does_not_move_it():
    before = Body(parse_file(BODY), seed=None).run(until=YEARS20)
    module = parse_file(BODY)
    assert _germ_layer_splits(module) == 16
    after = Body(module, seed=None).run(until=YEARS20)

    def census_of(body):
        out: dict[str, float] = defaultdict(float)
        for c in body.alive_at(YEARS20):
            out[c.cell_type] += c.count
        return out

    old, new = census_of(before), census_of(after)
    assert sum(new.values()) == pytest.approx(sum(old.values()), rel=1e-12)
    assert sum(old.values()) == pytest.approx(2.83e13, rel=5e-3)  # the count the program is solved to
    for kind in set(old) | set(new):
        assert new[kind] == pytest.approx(old[kind], rel=1e-9), kind
