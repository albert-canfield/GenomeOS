# SPDX-License-Identifier: AGPL-3.0-or-later
"""Stage 3's gate against its registration (v0.4 §9.5), and the ways it is allowed to fail.

The registered numbers are checked against their own closed forms first, so a typo in the registration
cannot pass for a prediction; then the committed programs against the registration; then the ten
clauses on the committed programs; then one break per clause that the verdict must notice.
"""

from __future__ import annotations

import copy
import importlib.util
from pathlib import Path

import pytest

from genomeos.lang import parse_file
from genomeos.runtime import metabolism_gate as reg

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "metabolism_gate_script", ROOT / "scripts" / "metabolism_gate.py"
)
assert spec and spec.loader
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


@pytest.fixture(scope="module")
def programs():
    return parse_file(ROOT / reg.PROGRAM), parse_file(ROOT / reg.RBC_PROGRAM)


def test_the_registered_numbers_agree_with_their_own_closed_forms():
    assert sum(c["rate_per_h"] * c["atp_per_event"] for c in reg.CONSUMERS.values()) == reg.DEMAND_PER_H
    assert (
        sum(c["rate_per_h"] * c["atp_per_event"] for c in reg.RBC_CONSUMERS.values()) == reg.RBC_DEMAND_PER_H
    )
    for arm, y in reg.PREDICTED_YIELD_AEROBIC.items():
        po = reg.ARMS[arm]
        assert 2 + 2 + 8 * po["po_nadh"] + 2 * po["po_fadh2"] + 2 * po["po_cytosolic"] == pytest.approx(y)
        assert reg.PREDICTED_THRESHOLD[arm] == pytest.approx(1 - 12 * reg.DEMAND_PER_H / (y * reg.CAPACITY))
        assert reg.PREDICTED_PASTEUR_RATIO[arm] == y / reg.PREDICTED_YIELD_ANAEROBIC
        avail = reg.ATP_STOCK + (reg.PREDICTED_MAX_SUPPLY[arm] - reg.DEMAND_PER_H) * reg.INTERVAL_H
        assert reg.PREDICTED_AVAILABLE[arm] == pytest.approx(avail)
    assert abs(reg.PREDICTED_RBC_YIELD - (2 - 2 * reg.BPG_SHUNT)) < 1e-12
    assert set(reg.PRIORITY_ORDER) == set(reg.CONSUMERS) and set(reg.MACROMOLECULE) <= set(reg.CONSUMERS)
    assert len(reg.PRE_REGISTRATION["clauses"]) == 10


def test_the_committed_programs_declare_what_was_registered(programs):
    cell, rbc = programs
    assert gate.disagreements(cell, rbc) == []
    changed = copy.deepcopy(cell)
    changed.entities["Glucose"].regeneration = 3e9
    assert any("Glucose uptake" in p for p in gate.disagreements(changed, rbc))


def test_every_registered_clause_is_met_on_the_committed_programs(programs):
    cell, rbc = programs
    got = {
        "1": gate.clause_yields(cell),
        "2": gate.clause_budget(cell),
        "3": gate.clause_pasteur(cell),
        "4": gate.clause_pasteur_uncapped(cell),
        "5": gate.clause_glycolysis_only(rbc),
        "6": gate.clause_glucose_removed(cell, rbc),
        "7": gate.clause_hierarchy(cell),
        "8": gate.clause_threshold(cell),
        "9": gate.clause_checkpoint(cell),
        "10": gate.clause_refusals(),
    }
    assert {k: v["ok"] for k, v in got.items()} == dict.fromkeys(got, True)


def test_a_scalar_policy_cannot_show_the_hierarchy(programs):
    """The registered negative: written as proportional, the cell loses the order clause 7 asks for."""
    cell, _ = programs
    assert not gate.clause_hierarchy(gate.with_policy(cell, "proportional"))["ok"]


def test_without_the_reserve_the_threshold_clause_fails(programs):
    """Capacity cut to what the demand needs: OXPHOS falls at the first heteroplasmy step, so the
    registered threshold is not where the run puts it."""
    cell, _ = programs
    tight = gate.with_param(cell, "oxphos_capacity_per_copy", 3.75e6)
    assert not gate.clause_threshold(tight)["ok"]


def test_a_glycolysis_that_read_oxygen_would_be_caught(programs, monkeypatch):
    _, rbc = programs
    real = gate.Metabolism.solve

    def leaky(self, arm):
        flux = real(self, arm)
        flux.supply *= 1.0 if self.o2_supply >= 3e9 else 0.99  # oxygen leaking into a glycolysis-only cell
        return flux

    monkeypatch.setattr(gate.Metabolism, "solve", leaky)
    assert not gate.clause_glycolysis_only(rbc)["ok"]
