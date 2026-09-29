# SPDX-License-Identifier: AGPL-3.0-or-later
"""The placement census (lane-census, 2026-09-29): the set-valued observation rule, the verdicts, the two
interval rules and the tables, on synthetic intervals. The run reads the benchmark and the compiled
programs; these tests do not."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

from genomeos.attribution import ablation as ab
from genomeos.attribution import measured as ms

ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location("placement_census", ROOT / "scripts/placement_census.py")
pc = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = pc
_SPEC.loader.exec_module(pc)
pa = pc.pa


def _pair(
    start, end, gene="G", cell="K562", split="training", effect=-0.3, significant=True, power=0.9, ds="D"
):
    return ms.CrispriPair(
        chrom="chr1",
        start=start,
        end=end,
        gene=gene,
        cell=cell,
        dataset=ds,
        reference="R",
        regulated=significant and effect < 0,
        significant=significant,
        effect_size=effect,
        p_adjusted=0.01 if significant else 0.9,
        split=split,
        power_at_effect_size_20=power,
    )


def _world(elements, registry_only=()):
    """elements: (start, end, gene, action) compiled predicted links on chr1, named E<i>; registry_only:
    (start, end) registry elements no link carries, named R<i>."""
    links = [pc.PredictedLink(i, "chr1", s, e, f"E{i}", g, a) for i, (s, e, g, a) in enumerate(elements)]
    reg = [(x.chrom, x.start, x.end, x.element) for x in links]
    reg += [("chr1", s, e, f"R{i}") for i, (s, e) in enumerate(registry_only)]
    return links, pa.ElementSet((x.chrom, x.start, x.end, x.index) for x in links), pa.ElementSet(reg)


def _observe(pairs, world, rule=pc.POLICY):
    links, link_set, registry = world
    return pc.observe(pairs, links, link_set, registry, pc.RULES[rule])


# --- the set-valued rule: one observation id never counts twice --------------------------------
def test_one_observation_over_two_links_is_one_supported_observation():
    # two 200 bp elements inside one 600 bp tested interval, both predicted to activate G
    world = _world([(100, 300, "G", "activates"), (300, 500, "G", "activates")])
    obs = _observe([_pair(0, 600)], world)
    (o,) = obs
    assert o.attached == (0, 1)
    assert o.candidates == ("E0", "E1")
    assert o.resolution == pc.AMBIGUOUS_SEVERAL
    assert o.verdict == pc.SUPPORTED
    table = pc.observation_table(obs)
    assert table["summary"]["supported"] == 1
    assert table["observations"] == 1
    assert table["by_verdict_and_resolution"][pc.SUPPORTED] == {pc.AMBIGUOUS_SEVERAL: 1}
    # the element level says two links, and that both rest on the one observation
    links, _, _ = world
    by_class, _, on, _ = pc.link_census(links, obs)
    assert by_class[(ab.SURVIVES, None)] == 2
    att = pc.element_attachments(links, obs, on)
    assert att["observation_link_pairs"] == 2
    assert att["supported_links"] == 2
    assert att["distinct_observations_behind_supported_links"] == 1
    assert att["supported_links_resting_only_on_ambiguous_observations"] == 2
    assert att["supported_links_with_a_unique_element_observation"] == 0


def test_a_repeated_observation_id_is_refused():
    world = _world([(100, 300, "G", "activates")])
    with pytest.raises(ValueError, match="twice"):
        _observe([_pair(100, 300), _pair(100, 300)], world)
    (o,) = _observe([_pair(100, 300)], world)
    with pytest.raises(ValueError, match="repeats"):
        pc.observation_table([o, o])


def test_every_observation_enters_each_table_once():
    world = _world(
        [(100, 300, "G", "activates"), (300, 500, "G", "activates"), (1000, 1200, "H", "inhibits")]
    )
    pairs = [
        _pair(0, 600),
        _pair(100, 300, cell="WTC11"),
        _pair(1000, 1200, gene="H", effect=0.3),
        _pair(1000, 1200, gene="G", significant=False, power=0.95),
        _pair(5000, 5500, significant=False, power=0.1),
    ]
    obs = _observe(pairs, world)
    t = pc.observation_table(obs)
    assert sum(t["by_verdict"].values()) == len(pairs)
    assert sum(t["by_resolution"].values()) == len(pairs)
    assert sum(sum(v.values()) for v in t["by_outcome_and_verdict"].values()) == len(pairs)
    s = t["summary"]
    assert sum(s.values()) == len(pairs)
    for kind, rows in pc.strata_tables(obs).items():
        assert sum(r["observations"] for r in rows.values()) == len(pairs), kind


# --- unique against ambiguous ------------------------------------------------------------------
def test_unique_needs_one_candidate_element():
    world = _world([(100, 300, "G", "activates")], registry_only=[(300, 500)])
    (alone,) = _observe([_pair(100, 300)], world)
    assert alone.resolution == pc.UNIQUE
    # the same link, but the tested interval also carries a registry element the link does not name
    (shared,) = _observe([_pair(0, 600)], world)
    assert shared.attached == (0,)
    assert shared.candidates == ("E0", "R0")
    assert shared.resolution == pc.AMBIGUOUS_ONE


def test_a_link_to_another_gene_is_a_candidate_but_not_attached():
    world = _world([(100, 300, "H", "activates")])
    (o,) = _observe([_pair(100, 300, gene="G")], world)
    assert o.candidates == ("E0",)
    assert o.compiled == (0,)
    assert o.attached == ()
    assert o.resolution == pc.UNATTACHED
    assert o.verdict == pc.UNATTACHED


# --- the verdicts ------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("actions", "kw", "want"),
    [
        (["activates"], {"effect": -0.3}, pc.SUPPORTED),
        (["activates"], {"effect": 0.3}, pc.REFUTED_OPPOSITE),
        (["inhibits"], {"effect": 0.3}, pc.SUPPORTED),
        (["inhibits"], {"effect": -0.3}, pc.REFUTED_OPPOSITE),
        (["activates", "inhibits"], {"effect": 0.3}, pc.SPLIT),
        (["activates"], {"significant": False, "power": 0.95}, pc.REFUTED_NULL),
        (["activates", "inhibits"], {"significant": False, "power": 0.95}, pc.REFUTED_NULL),
        (["activates"], {"significant": False, "power": 0.2}, pc.UNASSESSED_UNDERPOWERED),
        (["activates"], {"effect": float("nan")}, pc.UNASSESSED_MISSING),
    ],
)
def test_verdict(actions, kw, want):
    elements = [(100 + 200 * i, 300 + 200 * i, "G", a) for i, a in enumerate(actions)]
    (o,) = _observe([_pair(100, 100 + 200 * len(actions), **kw)], _world(elements))
    assert len(o.attached) == len(actions)
    assert o.verdict == want


# --- the two rules -----------------------------------------------------------------------------
def test_the_current_rule_attaches_what_crispri_index_attaches():
    elements = [
        (100, 300, "G", "activates"),
        (300, 500, "G", "activates"),
        (2000, 2250, "G", "inhibits"),
        (4000, 4200, "H", "activates"),
        (9000, 9300, "G", "activates"),
    ]
    pairs = [
        _pair(0, 600),  # 600 bp over two 200 bp elements: the reciprocal rule meets neither
        _pair(120, 320),
        _pair(1950, 2300, effect=0.2),
        _pair(4000, 4200, significant=False, power=0.9),
        _pair(4050, 4250, gene="G"),
        _pair(8800, 9700, significant=False, power=0.1),
    ]
    links, link_set, registry = world = _world(elements)
    obs = _observe(pairs, world, pc.CURRENT)
    by_class, by_detail, on, per_link = pc.link_census(links, obs)
    index = ab.CrispriIndex(pairs)
    direct = [ab.link_class(index.of(x.chrom, x.start, x.end), x.gene, x.action) for x in links]
    assert per_link == [d[:2] for d in direct]
    assert sum(by_class.values()) == len(links)
    for x in links:
        assert {pc.observation_id(p) for p in index.of(x.chrom, x.start, x.end)} == {
            o.id for o in on.get(x.index, ())
        }


def test_the_policy_keeps_every_current_attachment_and_the_difference_counts_wide_intervals():
    elements = [(100, 300, "G", "activates"), (300, 500, "G", "activates"), (3000, 3200, "G", "activates")]
    pairs = [
        _pair(100, 400),  # attached under both rules: to E0, and under the policy to E1 too
        _pair(0, 600),  # 600 bp: newly attached, to two links
        _pair(2600, 3600),  # 1,000 bp: newly attached, to one link
        _pair(7000, 7500),  # attached under neither
    ]
    links, _, _ = world = _world(elements)
    cur, pol = _observe(pairs, world, pc.CURRENT), _observe(pairs, world, pc.POLICY)
    lc = pc.link_census(links, cur)[3]
    lp = pc.link_census(links, pol)[3]
    d = pc.difference(cur, pol, lc, lp)
    assert d["attachments_lost_under_the_policy"] == {"candidate_sets": 0, "attached_links": 0}
    new = d["newly_attached_observations"]
    assert new["n"] == 2
    assert new["unique"] == 1
    assert new["ambiguous"] == 1
    assert new["by_width"] == {"gt_700": 1, "le_700": 1}
    edges = d["new_observation_link_attachments"]
    assert edges["n"] == 4  # two links for the 600 bp interval, one for the 1,000 bp, and E1 for 100-400
    assert edges["from_tested_intervals_over_700_bp"] == 1
    assert edges["on_newly_attached_observations"] == 3
    assert d["attached_under_both_rules"]["n"] == 1
    assert d["attached_under_both_rules"]["resolution_transitions"] == {
        f"{pc.UNIQUE}->{pc.AMBIGUOUS_SEVERAL}": 1
    }
    # three supported observations over five observation-link pairs: the count is of observations
    assert pc.observation_table(pol)["summary"]["supported"] == 3
    att = pc.element_attachments(links, pol, pc.link_census(links, pol)[2])
    assert att["observation_link_pairs"] == 5
    assert att["supported_links"] == 3
    assert att["distinct_observations_behind_supported_links"] == 3


def test_policy_rule_is_the_registered_one():
    assert pc.RULES[pc.CURRENT] is pa.current_rule
    assert pc.RULES[pc.POLICY] is pa.policy_rule
    assert pa.POLICY_THRESHOLD == ms.RECIPROCAL_OVERLAP
    assert not pa.current_rule(100, 300, 0, 600)
    assert pa.policy_rule(100, 300, 0, 600)


# --- the strata --------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("width", "coarse", "fine"),
    [
        (0, "zero_width", "zero_width"),
        (-5, "zero_width", "zero_width"),
        (52, "le_700", "1_to_74"),
        (75, "le_700", "75_to_350"),
        (500, "le_700", "351_to_500"),
        (700, "le_700", "501_to_700"),
        (701, "gt_700", "701_to_1000"),
        (2000, "gt_700", "1001_to_2000"),
        (4181, "gt_700", "over_2000"),
    ],
)
def test_width_classes(width, coarse, fine):
    assert pc.width_class(width) == coarse
    assert pc.width_fine(width) == fine


def test_a_zero_width_row_attaches_under_neither_rule():
    world = _world([(100, 300, "G", "activates")])
    for rule in pc.RULES:
        (o,) = _observe([_pair(200, 200)], world, rule)
        assert o.candidates == ()
        assert o.verdict == pc.UNATTACHED


def test_registration_holds_both_rules_the_set_valued_rule_and_the_stop():
    reg = pc.registration()
    assert set(reg["rules"]) == {pc.CURRENT, pc.POLICY}
    assert "never several validations" in reg["set_valued_rule"]
    assert "440,377" in reg["denominator"]
    assert any("recommendation" in x for x in reg["out_of_scope"])
    commits = {x["commit"] for x in pc.PRIOR_EXPOSURE}
    assert {"b7e4bf0", "dd8c49c", "6c0d39f", "0183ef5", "8f87bf6"} <= commits
