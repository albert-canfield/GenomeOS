# SPDX-License-Identifier: Apache-2.0
"""Stage 4: a division that partitions contents, and the checkpoint that guards it (v0.4 §8, §9.3).

The tests that matter most here are the ones that pin a *negative*: that partitioning is gene-blind,
so it cannot separate a stable protein from a short-lived one; that the default is the old behaviour,
so no committed program moves; and that a protein with no declared half-life is never given one.
"""

from __future__ import annotations

import math
import random

import pytest

from genomeos.ir import UNKNOWN, Event
from genomeos.lang import BioLangError, parse
from genomeos.runtime.division import Division, DivisionError

PROGRAM = """
module t.division
regime counted { units: copies; threshold: 50 }
protein Stable { half_life: 240 }
protein Short { half_life: 6 }
protein Forever { }
event divide { partition: contents; when: size >= 2 }
"""


def build(text: str = PROGRAM) -> Division:
    return Division.from_module(parse(text))


def test_default_partition_is_the_behaviour_every_program_already_had():
    """Stage 4 must be opt-in: a `divide` event that says nothing still duplicates."""
    assert Event(id="divide").partition == "duplicate"
    d = Division.from_module(parse("module t.d\nevent divide { rate: 1 /day }\n"))
    assert d.mode == "duplicate"
    got = d.split({"X": 100.0})
    assert got.a == got.b == {"X": 100.0}


def test_an_unknown_partition_mode_is_refused_by_name():
    with pytest.raises(BioLangError) as e:
        parse("module t.d\nevent divide { partition: sideways }\n")
    assert "partition 'sideways'" in str(e.value)
    assert "duplicate" in str(e.value) and "binomial" in str(e.value)


def test_an_unknown_key_on_an_event_is_still_refused():
    """The 2026-09-27 rule holds for the key stage 4 adds beside it."""
    with pytest.raises(BioLangError) as e:
        parse("module t.d\nevent divide { partitions: contents }\n")
    assert "no key 'partitions'" in str(e.value)


def test_partitioning_conserves_every_species_exactly():
    d = build()
    rng = random.Random(7)
    amounts = {"a": 3.0, "b": 17.0, "c": 1e6, "d": 1.0, "e": 0.0}
    for _ in range(200):
        got = d.split(amounts, rng)
        assert got.conserves(amounts)
    assert d.split({"c": 1e6}).rule == {"c": "halves"}
    assert d.split({"b": 17.0}).rule == {"b": "binomial"}


def test_binomial_refuses_half_a_molecule():
    d = build(PROGRAM.replace("partition: contents", "partition: binomial"))
    with pytest.raises(DivisionError, match="whole molecules"):
        d.split({"X": 17.5})


def test_a_protein_with_no_half_life_is_not_given_one():
    d = build()
    assert d.dilution_only == ("Forever",)
    assert d.half_lives["Forever"] is UNKNOWN
    assert d.decay({"Forever": 1e6}, 10_000.0) == {"Forever": 1e6}


def test_dilution_is_the_only_removal_that_reaches_it():
    """The positive clause of the gate: eight divisions take exactly 2^-8 of it, and nothing else."""
    d = build()
    got = d.run({"Forever": 1e6}, 8, 24.0, {"size": "2"})
    assert got["final"]["Forever"] == 1e6 * 2.0**-8
    duplicating = build(PROGRAM.replace("partition: contents", "partition: duplicate"))
    assert duplicating.run({"Forever": 1e6}, 8, 24.0, {"size": "2"})["final"]["Forever"] == 1e6


def test_dilution_cannot_separate_a_stable_protein_from_a_short_lived_one():
    """Registered before the instrument existed (runtime/division_gate.py) and pinned here."""
    start = {"Stable": 1e6, "Short": 1e6}
    partitioned = build().run(start, 8, 24.0, {"size": "2"})["final"]
    duplicated = build(PROGRAM.replace("partition: contents", "partition: duplicate")).run(
        start, 8, 24.0, {"size": "2"}
    )["final"]
    separation = math.log10(partitioned["Stable"] / partitioned["Short"])
    assert separation == math.log10(duplicated["Stable"] / duplicated["Short"])
    assert separation == pytest.approx(8 * 24 * math.log10(2) * (1 / 6 - 1 / 240), rel=1e-12)


def test_division_compresses_the_separation_in_a_steady_state():
    d = build()
    static = d.steady_state(1.0, 240.0, None) / d.steady_state(1.0, 6.0, None)
    dividing = d.steady_state(1.0, 240.0, 24.0) / d.steady_state(1.0, 6.0, 24.0)
    assert static == pytest.approx(40.0)
    assert static / dividing == pytest.approx(8.8)


def test_equal_half_lives_leave_the_two_proteins_identical():
    """The falsifier: any difference would be bookkeeping rather than turnover."""
    d = build(PROGRAM.replace("half_life: 6", "half_life: 240"))
    got = d.run({"Stable": 1e6, "Short": 1e6}, 8, 24.0, {"size": "2"})["final"]
    assert got["Stable"] == got["Short"]


def test_the_checkpoint_stops_the_division_and_not_the_clock():
    d = build()
    got = d.run({"Stable": 1e6}, 8, 24.0, {"size": "1"})
    assert got["divisions_fired"] == 0
    assert got["final"]["Stable"] == 1e6 * 2.0 ** (-8 * 24 / 240)
    assert "checkpoint unmet" in got["trace"][1]["why"]


def test_the_checkpoint_is_the_matcher_every_other_guard_uses():
    """Not a private mechanism: `a|b`, `absent` and `unknown` behave here as they do on a rule."""
    assert build().ready({"size": "3"})[0]
    d = build(PROGRAM.replace("when: size >= 2", "when: phase = S|G2"))
    assert d.ready({"phase": "G2"})[0] and not d.ready({"phase": "G1"})[0]
    d = build(PROGRAM.replace("when: size >= 2", "when: arrest = absent"))
    assert d.ready({})[0] and not d.ready({"arrest": "yes"})[0]
    d = build(PROGRAM.replace("when: size >= 2", "when: phase = unknown"))
    assert not d.ready({"phase": "S"})[0]


def test_a_species_with_no_removal_at_all_has_no_steady_state():
    with pytest.raises(DivisionError, match="no removal"):
        build().steady_state(1.0, UNKNOWN, None)


def test_partition_survives_the_ir_round_trip():
    """A key that did not round-trip would be a key that vanished on the way to a saved program."""
    from genomeos.ir import Module

    m = Module.from_dict(parse(PROGRAM).to_dict())
    event = next(e for e in m.events if e.id == "divide")
    assert event.partition == "contents"
    assert event.when == {"size": ">=2"}


def test_a_module_with_no_event_says_so_rather_than_dividing_nothing():
    with pytest.raises(DivisionError, match="declares no event"):
        Division.from_module(parse("module t.d\nprotein P { half_life: 1 }\n"))
