# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""BioLang v0.4 §7.5: a cell looks again when the quantity arrives, not when a clock ticks.

Until now the only way to be re-read was `cell_network`, which exists to step per-cell networks and
re-decides every cell as a side effect. The falsifier stated in §7.5 before this was built: if
re-decision is first-class, stepping a network must stop changing fates as a side effect, and a
program that names no integrated read must be identical under `crossings` and `none`.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from genomeos.lang import parse, parse_file
from genomeos.lang.parser import BioLangError
from genomeos.runtime.body import Body

ROOT = Path(__file__).resolve().parent.parent
WORM = ROOT / "data" / "organisms" / "celegans" / "embryo_factors.bio"
PLAIN_WORM = ROOT / "data" / "organisms" / "celegans" / "embryo.bio"

WAITING = """
module toy
import bio.std.development
cell_type Ready { parent: PostMitotic }
organism T { root: P0; cell_type: Blastomere; factors: M; observe: fates }
stage S { from: 0 min }
decision wait { action: quiesce; when: cell = P0, M.exposure(cell) = <25 }
decision ready { action: differentiate; when: cell = P0, M.exposure(cell) >= 25; to: Ready }
"""


def _when_decided(source: str, recheck: str, until: float = 100.0) -> tuple[list[float], Body]:
    """The times at which any cell decided, so the crossing can be checked to the minute."""
    import genomeos.runtime.body as body_module

    times: list[float] = []
    original = body_module.Body._decide

    def watched(self, c, born=True):
        times.append(self.time)
        return original(self, c, born)

    body_module.Body._decide = watched
    try:
        body = Body(parse(source), seed=None, recheck=recheck).run(until=until)
    finally:
        body_module.Body._decide = original
    return times, body


def test_a_cell_wakes_at_the_instant_its_read_reaches_the_threshold():
    times, body = _when_decided(WAITING, "crossings")
    assert body.cells["P0"].cell_type == "Ready"
    assert times[-1] == pytest.approx(25.0)  # exactly when exposure reaches 25, not a tick before or after
    assert body.summary()["rechecks"] == 1


def test_without_crossings_the_cell_never_looks_again():
    _, body = _when_decided(WAITING, "none")
    assert body.cells["P0"].cell_type == "Blastomere"
    assert body.summary()["rechecks"] == 0


def test_a_mean_over_a_window_of_no_length_is_the_value_in_force():
    """At the instant a cell is born its mean is 0/0; the limit is the rate now in force, so a cell
    reading `mean(cell)` at birth is taking the instantaneous read and needs no wake-up for it."""
    source = WAITING.replace("M.exposure(cell) = <25", "M.mean(cell) = <0.5").replace(
        "M.exposure(cell) >= 25", "M.mean(cell) >= 0.5"
    )
    _, body = _when_decided(source, "crossings")
    assert body.cells["P0"].cell_type == "Ready"
    assert body.summary()["rechecks"] == 0


def test_the_regime_declares_it_and_refuses_anything_else():
    assert Body(parse(WAITING + "regime r { recheck: none }\n"), seed=None).recheck == "none"
    assert Body(parse(WAITING), seed=None).recheck == "crossings"  # the default
    with pytest.raises(BioLangError, match="recheck must be one of"):
        parse(WAITING + "regime r { recheck: sometimes }\n")


def test_a_program_that_names_no_read_schedules_nothing():
    plain = parse_file(PLAIN_WORM)
    body = Body(plain, seed=None, recheck="crossings").run(until=200)
    assert body._terms == {} and body.summary()["rechecks"] == 0


def _terminal_fates(cadence: int, recheck: str) -> tuple[dict[str, tuple[str, str]], dict]:
    module = parse_file(WORM)
    module.organism.cell_network = cadence
    body = Body(module, seed=None, recheck=recheck).run(until=800)
    return {n: (c.cell_type, c.terminal_name) for n, c in body.cells.items() if c.terminal}, body.summary()


def test_stepping_a_network_no_longer_changes_which_fates_are_taken():
    """The falsifier of §7.5, stated before it was built: the worm's terminal fates and their names
    must be the same whether or not a network is being stepped."""
    still, still_summary = _terminal_fates(0, "crossings")
    stepped, stepped_summary = _terminal_fates(6, "crossings")
    assert stepped == still
    assert still_summary["rechecks"] > 0 and stepped_summary["network_steps"] > 0


def test_a_cell_re_read_at_a_crossing_keeps_its_name():
    """A decision that lands on the type the cell already holds revises nothing, so it must not take
    the cell's terminal name away — the §7.4 anchor-cell bug, reachable again from the other side."""
    named = {}
    for recheck in ("none", "crossings"):
        body = Body(parse_file(WORM), seed=None, recheck=recheck).run(until=800)
        named[recheck] = {n: c.terminal_name for n, c in body.cells.items() if c.terminal_name}
    assert named["crossings"] == named["none"]
    assert len(named["none"]) > 500
