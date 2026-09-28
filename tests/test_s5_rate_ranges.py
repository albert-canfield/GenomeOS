# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item 12 S5 (registered in 0ea0a90): simulation on transferred rates is labelled on its output.

Negatives first: a rate whose species nobody stated is labelled unstated and never assumed; a
compiled program run on measured rates without a half-life names the defaulted half-life, because
its level is in mixed units; a program with no transferred rate carries no transfer label. Then the
label on every path (the bridge's result, the gene's attributes, the trajectory, the audit's counts),
the combination rule's irrelevance to a single-mechanism gene, and the propagation's own arithmetic:
the switch band's interior, an inhibition that can never be switched off, the rank count, the gate.
"""

from __future__ import annotations

import math

import pytest

from genomeos.attribution import bridge
from genomeos.lang.parser import parse
from genomeos.runtime import grn
from genomeos.runtime.grn import NetworkRuntime

K562 = {"cell_type": "K562"}
T_MEASURED = 12.0
MOUSE = {"source": "schwanhausser2011", "species_cell": "Mus musculus, NIH 3T3 fibroblast"}


def _program(lfc: float = -1.0, half_life: float | None = 4.0) -> object:
    param = (
        f'param mrna_half_life = {half_life} h {{ evidence: experimental "Schwanhausser 2011, mouse" }}\n'
        if half_life is not None
        else ""
    )
    return parse(
        param + "element E1 { class: enhancer }\n"
        "gene GENE1 { symbol: GENE1 }\n"
        f"rule E1 {'activates' if lfc < 0 else 'inhibits'} GENE1 {{ strength: {abs(lfc)};"
        f' when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" effect {lfc:+}'
        " log2 fold change on deletion }\n",
        "t",
    )


def _run(module, clamp=None):
    return NetworkRuntime(module, context=K562, strict=True).run(
        hours=40.0, dt=0.01, clamp=clamp or {"E1": 1.0}, record_every=1000
    )


# ---- negatives -------------------------------------------------------------------------------------


def test_a_rate_whose_species_nobody_stated_is_labelled_unstated_not_assumed():
    b = bridge.parameterize(_program(), K562, rates={"GENE1": T_MEASURED})
    lab = b.transferred["GENE1"]
    assert lab["species_cell"] == "unstated" and lab["source"] == "unstated"
    tr = _run(b.module).transferred
    assert [(t.subject, t.species_cell) for t in tr] == [("GENE1", "unstated")]
    assert not any("Homo" in t.species_cell or "Mus" in t.species_cell for t in tr)


def test_a_compiled_program_on_measured_rates_names_its_defaulted_half_life():
    """S5_RUNTIME_GAP: T in molecules over an a.u. decay is mixed units, and the trajectory says so."""
    b = bridge.parameterize(
        _program(half_life=None), K562, rates={"GENE1": T_MEASURED}, rate_provenance={"GENE1": MOUSE}
    )
    traj = _run(b.module)
    assert traj.defaulted == ["mrna_half_life"]
    assert traj.transferred and traj.provenance()["defaulted"] == ["mrna_half_life"]


def test_a_program_without_a_transferred_rate_carries_no_transfer_label():
    b = bridge.parameterize(_program(), K562, genes={"GENE1": {"basal_rate": 1.0, "max_rate": 10.0}})
    assert b.transferred == {}
    traj = _run(b.module)
    assert traj.transferred == [] and traj.defaulted == []


def test_an_unfitted_gene_gets_no_label():
    """A gene reported unresolved is not a simulation output, so nothing is labelled on its behalf."""
    b = bridge.parameterize(
        _program(lfc=1.5), K562, rates={"GENE1": T_MEASURED}, rate_provenance={"GENE1": MOUSE}
    )
    assert [i.reason for i in b.unresolved] == ["rate_split_infeasible"] and b.transferred == {}


# ---- the label on every path -----------------------------------------------------------------------


def test_the_label_travels_from_the_bridge_to_the_gene_to_the_trajectory():
    b = bridge.parameterize(_program(), K562, rates={"GENE1": T_MEASURED}, rate_provenance={"GENE1": MOUSE})
    lab = b.transferred["GENE1"]
    assert lab == {
        "quantity": grn.TRANSFERRED_RATE_QUANTITY,
        "tier": "measured",
        "source": "schwanhausser2011",
        "species_cell": "Mus musculus, NIH 3T3 fibroblast",
        "used_in": "K562",
        "status": bridge.SIMULATION_STATUS,
    }
    g = b.module.entities["GENE1"]
    assert [g.attrs[k] for k in grn.RATE_PROVENANCE_ATTRS] == [
        "schwanhausser2011",
        "Mus musculus, NIH 3T3 fibroblast",
        "measured",
    ]
    for clamp in ({"E1": 1.0}, {"E1": 0.0}):
        (t,) = _run(b.module, clamp).transferred
        assert (t.subject, t.source, t.species_cell, t.tier, t.used_in) == (
            "GENE1",
            "schwanhausser2011",
            "Mus musculus, NIH 3T3 fibroblast",
            "measured",
            "K562",
        )


def test_the_audit_mode_labels_without_building():
    b = bridge.parameterize(
        _program(), K562, build=False, rates={"GENE1": T_MEASURED}, rate_provenance={"GENE1": MOUSE}
    )
    assert b.transferred["GENE1"]["species_cell"].startswith("Mus musculus")


def test_a_borrowed_median_is_labelled_as_borrowed_with_its_own_provenance():
    lent = {"source": "genome median of schwanhausser2011", "species_cell": MOUSE["species_cell"]}
    b = bridge.parameterize(_program(), K562, rates={}, borrowed_rate=1.76, borrowed_provenance=lent)
    assert b.transferred["GENE1"]["tier"] == "borrowed_median"
    assert b.transferred["GENE1"]["source"].startswith("genome median")
    (t,) = _run(b.module).transferred
    assert t.tier == "borrowed_median"


def test_the_audit_counts_every_simulable_pair_as_labelled():
    from scripts.bridge_audit import audit_module  # noqa: PLC0415 - script, not a package

    out = audit_module(
        _program(),
        rates={"GENE1": T_MEASURED},
        borrowed=1.0,
        species={"GENE1": MOUSE["species_cell"]},
        provenance=({"GENE1": MOUSE}, {}),
    )
    for tier in ("measured", "borrowed"):
        c = out["by_tier"][tier]
        assert c["labelled_pairs"] == c["simulable_gene_cell_pairs"] == 1
        assert c["labelled_pairs_species_unstated"] == 0 and c["labelled_pairs_on_a_human_rate"] == 0


# ---- the combination rule and the model ------------------------------------------------------------


@pytest.mark.parametrize("rule", grn.ACTIVATOR_RULES)
def test_a_single_mechanism_gene_runs_the_same_under_every_activator_rule(rule, monkeypatch):
    b = bridge.parameterize(_program(), K562, rates={"GENE1": T_MEASURED}, rate_provenance={"GENE1": MOUSE})
    base = _run(b.module).final()["GENE1.mRNA"]
    monkeypatch.setattr(grn, "ACTIVATOR_COMBINATION", rule)
    assert _run(b.module).final()["GENE1.mRNA"] == pytest.approx(base, rel=1e-12)


def test_the_closed_form_is_the_runtime_for_a_one_gene_pair():
    from scripts.s5_rate_ranges import closed_form, runtime  # noqa: PLC0415

    for rho in (0.5, 0.9, 1.4):
        got = runtime(rho, 3.0, 2.5, "K562", MOUSE)
        want = closed_form(rho, 3.0, 2.5)
        assert got["intact"] == pytest.approx(want["intact"], rel=0.01)
        assert got["removed"] == pytest.approx(want["removed"], rel=0.01)
        assert got["response_time"] == pytest.approx(want["response_time"], rel=0.02)
        assert all(lab["transferred"][0]["used_in"] == "K562" for lab in got["labels"])


# ---- the propagation's arithmetic ------------------------------------------------------------------


def test_the_switch_band_interior_is_seen_even_when_every_corner_says_no():
    """S5_AMENDMENT (1): corners at 0.5 and 10 molecules both say 'no'; the band 1-4 lies between."""
    from scripts.s5_rate_ranges import classify, switch_stable  # noqa: PLC0415

    rho = 0.25
    vals = [{"intact": x, "removed": rho * x, "response_time": 1.0} for x in (0.5, 10.0)]
    assert not any(v["intact"] >= 1 and v["removed"] < 1 for v in vals)
    assert classify(rho, vals)["switch"] is False
    assert switch_stable(rho, 1.5, 3.0) and switch_stable(rho, 5.0, 9.0) and switch_stable(rho, 0.1, 0.9)


def test_an_inhibition_can_never_be_switched_off():
    """S5_AMENDMENT (2): removal raises an inhibited gene, so its switch answer is always stable."""
    from scripts.s5_rate_ranges import classify  # noqa: PLC0415

    rho = 1.5
    vals = [{"intact": x, "removed": rho * x, "response_time": 1.0} for x in (0.2, 30.0)]
    c = classify(rho, vals)
    assert c["on_off"] is False and c["switch"] is True


def test_a_level_is_stable_only_inside_two_fold_and_direction_and_fold_are_the_input():
    from scripts.s5_rate_ranges import classify, closed_form  # noqa: PLC0415

    narrow = [closed_form(0.5, x, 2.0) for x in (1.0, 1.9)]
    wide = [closed_form(0.5, x, 2.0) for x in (1.0, 2.1)]
    assert classify(0.5, narrow)["level"] and not classify(0.5, wide)["level"]
    c = classify(0.5, wide)
    assert c["direction"] and c["fold"]


def _registered_ranges() -> dict:
    reg = bridge.S5_RANGES
    return {
        "w_T": reg["R1_measurement"]["w_T"],
        "w_t_half": reg["R1_measurement"]["w_t_half"],
        "q05": reg["R2_transfer_to_a_human_cell"]["q05"],
        "q50": reg["R2_transfer_to_a_human_cell"]["q50"],
        "q95": reg["R2_transfer_to_a_human_cell"]["q95"],
        "R3_T": reg["R3_no_transfer"]["T"],
        "R3_t_half": reg["R3_no_transfer"]["t_half"],
    }


def test_the_ranges_points_follow_the_registration():
    from scripts.s5_rate_ranges import points  # noqa: PLC0415

    pts = points(2.0, 10.0, _registered_ranges())
    assert [len(pts[k]) for k in ("R0", "R1", "R2", "R3")] == [1, 4, 4, 4]
    for name, rate, t in pts["R2"]:
        if "copy number" in name:
            assert rate * t == pytest.approx(20.0)  # the mouse copy number carries over: the level stays
        else:
            assert rate == 2.0 and t < 10.0


def test_the_gate_refuses_a_range_that_is_not_the_registered_one():
    from scripts.s5_rate_ranges import gate  # noqa: PLC0415

    rg = _registered_ranges()
    assert gate(rg) == []
    assert gate({**rg, "q95": rg["q95"] + 0.01})


def test_rank_comparisons_count_only_ranges_that_cannot_swap():
    from scripts.s5_rate_ranges import stable_comparisons  # noqa: PLC0415

    by_cell = {"K562": [(0.0, 1.0), (0.5, 1.5), (3.0, 4.0)], "HepG2": [(0.0, 1.0)]}
    assert stable_comparisons(by_cell) == (2, 3)  # the third lies above both; the first two overlap


def test_nothing_is_ever_called_validated():
    assert bridge.SIMULATION_STATUS == "assumption-dependent simulation, not measured human dynamics"
    assert "NOT validated" in bridge.VALIDATED_HUMAN_KINETICS
    assert math.isclose(2.0 ** bridge.S5_RANGES["R2_transfer_to_a_human_cell"]["q95"], 0.6392, abs_tol=1e-4)
