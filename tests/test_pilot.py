# SPDX-License-Identifier: AGPL-3.0-or-later
"""The coherence pilot (item 13 C1-C3, item 12 S3): the energy, the debugger, boundaries, families and
the next measurement, on toy neighbourhoods that are not the registered gate's instances."""

from __future__ import annotations

import math

import pytest

from genomeos.attribution import joint_pretest as jp
from genomeos.attribution import pilot as pl

C = "chrT"


def _priors(tss, links=None, k562=(), hepg2=(), resp=None):
    return pl.Priors(
        C, dict(tss), dict(links or {}), pl.Peaks({"K562": list(k562), "HepG2": list(hepg2)}), resp
    )


def _obs(rows):
    out = []
    for kind, s, e, pos, gene, ctx in rows:
        out.append(pl.make_obs(len(out), kind, C, s, e, pos, gene=gene, ctx=ctx, source="toy"))
    return out


def test_a_spanning_observation_is_explained_by_either_block_once():
    tss = {"GA": 90_000}
    obs = _obs([("crispri", 0, 2000, True, "GA", "K562")])
    pri = _priors(tss, k562=[(0, 2000)])
    h = pl.Hood(C, [("A", 100, 500), ("B", 1200, 1600)], obs, pri)
    st = h.s0()
    assert h.violated(st) == [0]
    a, _ = h.step(st, h.energy(st), ("target", "A", "GA"))
    b, _ = h.step(st, h.energy(st), ("target", "B", "GA"))
    assert h.fits(a) and h.fits(b)
    assert h.signature(a) == h.signature(b)


def test_step_energy_equals_a_full_recomputation():
    tss = {"GA": 50_000, "GB": 150_000}
    obs = _obs(
        [
            ("crispri", 0, 1000, True, "GB", "K562"),
            ("crispri", 0, 1000, False, "GA", "K562"),
            ("reporter", 0, 500, False, "", "K562"),
            ("gtex", 400, 900, True, "GB", ""),
        ]
    )
    h = pl.Hood(
        C,
        [("A", 0, 1000), ("B", 1200, 1500)],
        obs,
        _priors(tss, {"A": ("GA", 0.7, "K562")}, k562=[(0, 1500)]),
    )
    st = h.s0()
    e = h.energy(st)
    for m in [
        ("target", "A", "GB"),
        ("act", "A", "K562"),
        ("split", "A", 500),
        ("merge", "A", "B"),
        ("excuse", "A", "reporter"),
    ]:
        s2, e2 = h.step(st, e, m)
        assert e2 == pytest.approx(h.energy(s2)), m


def test_a_split_conserves_the_prior_and_costs_fragmentation():
    tss = {"GA": 500_000}
    h = pl.Hood(C, [("A", 0, 1000)], [], _priors(tss, {"A": ("GA", 0.5, "K562")}, k562=[(0, 1000)]))
    st = h.s0()
    e = h.energy(st)
    s2, e2 = h.step(st, e, ("split", "A", 500))
    # the distance term moves by the halves' midpoints only; the rest is the split cost
    assert e2 - e == pytest.approx(pl.W["split"], abs=0.01)


def test_a_merge_does_not_lend_a_compiled_target_to_the_neighbour():
    tss = {"GA": 60_000, "GB": 80_000}
    links = {"A": ("GA", 1.0, "K562"), "B": ("GB", 1.0, "K562")}
    h = pl.Hood(C, [("A", 0, 300), ("B", 400, 700)], [], _priors(tss, links, k562=[(0, 700)]))
    st = h.s0()
    e = h.energy(st)
    _, em = h.step(st, e, ("merge", "A", "B"))
    _, er = h.step(st, e, ("target", "B", "GA"))
    # merging B under A's label costs B's own compiled target, like retargeting B would, plus the merge
    assert em - e == pytest.approx(er - e + pl.W["merge"], abs=0.05)


def test_the_debugger_finds_a_repair_that_needs_two_simultaneous_changes():
    # a toy swap: each block alone cannot move without breaking a screen that spans both
    tss = {"GA": 30_000, "GB": 45_000}
    links = {"A": ("GA", 1.2, "K562"), "B": ("GB", 1.2, "K562")}
    obs = _obs(
        [
            ("crispri", 0, 1400, True, "GA", "K562"),
            ("crispri", 0, 1400, True, "GB", "K562"),
            ("crispri", 100, 500, False, "GA", "K562"),
            ("crispri", 100, 500, True, "GB", "K562"),
            ("crispri", 900, 1300, True, "GA", "K562"),
            ("crispri", 900, 1300, False, "GB", "K562"),
        ]
    )
    h = pl.Hood(C, [("A", 100, 500), ("B", 900, 1300)], obs, _priors(tss, links, k562=[(0, 1400)]))
    greedy = pl.search(h, pairs=False)
    assert greedy.best.labels["A"].target == "GA"  # stuck: no single change helps
    out = pl.search(h)
    assert out.best.labels["A"].target == "GB" and out.best.labels["B"].target == "GA"
    assert h.fits(out.best)
    assert out.status == "resolved"
    assert {c["kind"] for c in out.committed} == {"target"}


def test_a_split_is_proposed_only_where_evidence_changes():
    tss = {"GA": 20_000, "GB": 70_000}
    obs = _obs([("crispri", 0, 700, True, "GB", "K562"), ("crispri", 700, 1500, True, "GA", "K562")])
    h = pl.Hood(C, [("A", 0, 1500)], obs, _priors(tss, {"A": ("GA", 0.5, "K562")}, k562=[(-100, 1600)]))
    st = h.s0()
    splits = [m for m in h.moves(st, ["A"]) if m[0] == "split"]
    assert splits == [("split", "A", 700)]


def test_indistinguishable_alternatives_form_one_family_and_the_pilot_abstains():
    tss = {"GA": 40_000}
    obs = _obs([("crispri", 0, 1500, True, "GA", "K562")])
    h = pl.Hood(C, [("A", 100, 500), ("B", 900, 1300)], obs, _priors(tss, k562=[(0, 1500)]))
    out = pl.search(h)
    assert out.status == "abstained"
    fam = out.families[0]
    assert fam["size"] >= 2
    nm = out.next_measurement
    assert nm is not None and nm["measure"] == "crispri" and nm["gene"] == "GA"
    assert nm["interval"] in (f"{C}:100-500", f"{C}:900-1300")


def test_the_smallest_conflicting_set_names_the_block_boundary():
    tss = {"GA": 20_000, "GB": 70_000}
    obs = _obs(
        [
            ("crispri", 0, 700, True, "GB", "K562"),
            ("crispri", 700, 1500, True, "GA", "K562"),
            ("reporter", 0, 1500, True, "", "K562"),
        ]
    )
    h = pl.Hood(C, [("A", 0, 1500)], obs, _priors(tss, {"A": ("GA", 0.5, "K562")}, k562=[(-100, 1600)]))
    st = h.s0()
    (i,) = h.violated(st)
    m = pl.mus(h, st, i)
    assert m["minimal"]
    # the smallest set: the screen and the label it contradicts
    assert m["size"] == 2 and m["assumptions"] == ["A targets GA"]
    # the conflict among the observations themselves: two screens, one block, one target
    among = m["conflict_among_observations"]
    assert len(among["observations"]) == 2  # the reporter tile is not part of the conflict
    assert "a block has one target" in among["rules"]
    assert any("boundary" in a for a in among["assumptions"])
    assert not any("targets" in a for a in among["assumptions"])  # the target itself may change


def test_a_pair_score_is_a_contrast_with_no_target_never_a_share_over_genes():
    tss = {"GA": 30_000, "GB": 31_000}
    obs = _obs([("crispri", 0, 1000, True, "GA", "K562")])
    h1 = pl.Hood(C, [("A", 0, 1000)], obs, _priors({"GA": 30_000}, k562=[(0, 1000)]))
    h2 = pl.Hood(C, [("A", 0, 1000)], obs, _priors(tss, k562=[(0, 1000)]))
    o1, o2 = pl.search(h1), pl.search(h2)
    s1 = pl.root_score(h1, o1, "A", 0, 1000, gene="GA", ctx="K562")
    s2 = pl.root_score(h2, o2, "A", 0, 1000, gene="GA", ctx="K562")
    assert s1 == pytest.approx(s2)  # a second candidate gene does not dilute the first


def test_the_independent_mode_has_no_gene_coupling_and_prior_mode_no_observations():
    resp = pl.Responsiveness({"GA": [(50_000, 51_000)]}, {})
    obs = _obs([("crispri", 0, 1000, True, "GA", "K562")])
    pri = _priors({"GA": 30_000}, k562=[(0, 1000)], resp=resp)
    joint = pl.Hood(C, [("A", 0, 1000)], obs, pri)
    prior = pl.Hood(C, [("A", 0, 1000)], obs, pri, mode="prior")
    assert prior.obs == []
    no_resp = _priors({"GA": 30_000}, k562=[(0, 1000)])
    ind = pl.Hood(C, [("A", 0, 1000)], obs, no_resp, mode="independent")
    assert joint.tlogw("A", "GA") - ind.tlogw("A", "GA") == pytest.approx(pl.W["responsive"])


def test_crispri_is_never_excused_and_fixed_blocks_never_move():
    tss = {"GA": 30_000}
    obs = _obs([("crispri", 0, 1000, True, "GA", "K562"), ("reporter", 0, 1000, True, "", "HepG2")])
    h = pl.Hood(C, [("A", 0, 1000), ("F", 2000, 2500)], obs, _priors(tss), fixed=["F"])
    st = h.s0()
    moves = h.moves(st, ["A", "F"])
    assert ("excuse", "A", "reporter") in moves
    assert not any(m[0] == "excuse" and m[2] == "crispri" for m in moves)
    assert not any(m[1] == "F" for m in moves)


def test_r8s_two_transformations_are_retired_and_its_reading_kept_beside_the_correction():
    assert jp.RETIRED == ("element_competition", "gene_budget")
    assert "cannot beat per-block scoring" in jp.READINGS["neither_passes"]  # the original, kept
    assert "does not show" in jp.READING_CORRECTIONS["neither_passes"]


def test_the_weights_are_log_odds_and_the_registration_is_complete():
    reg = pl.registration()
    assert (
        reg["weights"]["crispri_decrease"]
        > reg["weights"]["crispri_null"]
        > reg["weights"]["gtex_association"]
    )
    assert set(reg) >= {"weights", "search", "contexts", "measurement_cost", "retired_never_used"}
    assert all(math.isfinite(v) for v in pl.W.values())
