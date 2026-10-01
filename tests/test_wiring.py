# SPDX-License-Identifier: AGPL-3.0-or-later
"""The wiring diagnostic's matching, on synthetic inputs only."""

from __future__ import annotations

import math
import random

import pytest

from genomeos.attribution import crispri, wiring


def pair(
    chrom: str,
    start: int,
    gene: str,
    distance: float,
    regulated: bool = False,
    covered: bool = True,
    answered: bool = True,
    drop: float = 0.1,
    top: float = 0.0,
    activity: float = 1.0,
) -> crispri.Pair:
    p = crispri.Pair(chrom, start, start + 500, gene, "K562", "synthetic", distance, 1.0, 1.0, regulated)
    p.covered = covered
    p.features = {
        "log_distance": math.log(max(crispri.MIN_DISTANCE, distance)),
        "log_activity": activity,
        "activity_over_distance": activity - math.log(max(crispri.MIN_DISTANCE, distance)),
        "top_target": top,
        "deletion_drop": drop if covered and answered else 0.0,
        "deletion_answered": float(covered and answered),
        "covered": float(covered),
    }
    return p


def synthetic(n_elements: int = 30, genes_per_element: int = 8, seed: int = 1) -> list[crispri.Pair]:
    """Elements measured against several genes at spread distances; some pairs outside the population."""
    rng = random.Random(seed)
    out = []
    for e in range(n_elements):
        chrom = f"chr{1 + e % 3}"
        start = 1_000_000 * (e + 1)
        act = rng.uniform(0, 3)
        for g in range(genes_per_element):
            d = 10 ** rng.uniform(3.2, 5.8)
            out.append(
                pair(
                    chrom,
                    start,
                    f"G{e}_{g}",
                    d,
                    regulated=(g == 0),
                    covered=(g != genes_per_element - 1),
                    answered=(g != genes_per_element - 2),
                    drop=rng.uniform(0, 1),
                    top=float(g == 0),
                    activity=act,
                )
            )
    return out


def expression_for(pairs: list[crispri.Pair]) -> dict[str, float]:
    return {p.gene: (i % 7) / 2 for i, p in enumerate(pairs) if i % 5}


# --- derangement -------------------------------------------------------------------------------------------


def test_derange_never_gives_a_link_its_own_label():
    rng = random.Random(0)
    for labels in (
        ["a", "b"],
        ["a", "b", "c"],
        ["a", "a", "b", "b"],
        ["a", "a", "b", "c", "d"],
        list("abcdefgh"),
    ):
        for _ in range(50):
            src = wiring.derange(labels, rng)
            assert sorted(src) == list(range(len(labels)))
            assert all(labels[src[k]] != labels[k] for k in range(len(labels)))


def test_derange_refuses_when_no_assignment_exists():
    rng = random.Random(0)
    assert wiring.derange(["a"], rng) is None
    assert wiring.derange([], rng) is None
    assert wiring.derange(["a", "a"], rng) is None
    assert wiring.derange(["a", "a", "b"], rng) is None  # "a" fills more than half


def test_derange_constructive_fallback_is_valid():
    rng = random.Random(3)
    labels = ["a", "a", "a", "b", "b", "c"]
    for _ in range(20):
        src = wiring.derange(labels, rng, tries=0)
        assert sorted(src) == list(range(len(labels)))
        assert all(labels[src[k]] != labels[k] for k in range(len(labels)))


# --- one rewiring --------------------------------------------------------------------------------


@pytest.mark.parametrize("null", list(wiring.NULLS))
def test_rewire_respects_the_matching(null):
    pairs = synthetic()
    links = wiring.links_of(pairs, expression_for(pairs))
    groups, distinct = wiring.keys_for(links, null, expression_classes=False)
    pos = {lk.index: k for k, lk in enumerate(links)}
    for d in range(10):
        rw = wiring.rewire(links, groups, distinct, 0.3, wiring.rng_for("t", null, "S", d))
        assert sorted(rw) == sorted(rw.values())  # a permutation of the moved links
        for t, s in rw.items():
            assert groups[pos[t]] == groups[pos[s]]
            assert distinct[pos[t]] != distinct[pos[s]]
            assert abs(links[pos[t]].position - links[pos[s]].position) < 0.3


def test_element_kept_attaches_the_same_elements_value_on_another_gene():
    pairs = synthetic()
    links = wiring.links_of(pairs, {})
    groups, distinct = wiring.keys_for(links, "element_kept", expression_classes=False)
    rw = wiring.rewire(links, groups, distinct, 0.5, random.Random(5))
    assert rw
    for t, s in rw.items():
        assert wiring.element_of(pairs[t]) == wiring.element_of(pairs[s])
        assert pairs[t].gene != pairs[s].gene


def test_rewiring_is_deterministic_and_never_reads_labels():
    pairs = synthetic()
    links = wiring.links_of(pairs, {})
    groups, distinct = wiring.keys_for(links, "element_kept", expression_classes=False)
    a = wiring.rewire(links, groups, distinct, 0.3, wiring.rng_for("t", "element_kept", "S1", 4))
    b = wiring.rewire(links, groups, distinct, 0.3, wiring.rng_for("t", "element_kept", "S1", 4))
    assert a == b
    for p in pairs:
        p.regulated = not p.regulated
    flipped = wiring.links_of(pairs, {})
    c = wiring.rewire(flipped, groups, distinct, 0.3, wiring.rng_for("t", "element_kept", "S1", 4))
    assert a == c


def test_population_is_covered_pairs_with_a_cached_value():
    pairs = synthetic(n_elements=2, genes_per_element=8)
    links = wiring.links_of(pairs, {})
    assert len(links) == 2 * 6  # one uncovered and one unanswered pair per element stay outside
    assert all(wiring.in_population(pairs[lk.index]) for lk in links)


# --- fair scoring: the evaluated pairs and their labels never change ---------------------------------------


def test_attach_changes_only_the_block_of_moved_pairs():
    pairs = synthetic()
    real_rows = [dict(p.features) for p in pairs]
    real_labels = [p.regulated for p in pairs]
    real_ids = [(p.chrom, p.start, p.end, p.gene) for p in pairs]
    links = wiring.links_of(pairs, {})
    for null in wiring.NULLS:
        groups, distinct = wiring.keys_for(links, null, expression_classes=False)
        for d in range(5):
            rw = wiring.rewire(links, groups, distinct, 0.3, wiring.rng_for("t", null, "S", d))
            rows = wiring.attach(pairs, rw)
            # the same pairs, in the same order, with the same labels, in real and every rewiring
            assert len(rows) == len(pairs)
            assert [(p.chrom, p.start, p.end, p.gene) for p in pairs] == real_ids
            assert [p.regulated for p in pairs] == real_labels
            for i, (row, real) in enumerate(zip(rows, real_rows, strict=True)):
                for col in real:
                    if col in wiring.BLOCK and i in rw:
                        assert row[col] == real_rows[rw[i]][col]
                    else:
                        assert row[col] == real[col]
            # pairs outside the population and immovable links keep their real values
            for i, p in enumerate(pairs):
                if not wiring.in_population(p) or i not in rw:
                    assert rows[i] == real_rows[i]
    assert [dict(p.features) for p in pairs] == real_rows  # attach never writes into the pairs


# --- diagnostics --------------------------------------------------------------------------------


def test_alternatives_are_the_admissible_partners():
    links = [
        wiring.Link(0, ("c", 1, 2), "A", "c", 4.00, 1.0, None, True),
        wiring.Link(1, ("c", 1, 2), "B", "c", 4.15, 1.0, None, False),
        wiring.Link(2, ("c", 1, 2), "C", "c", 4.50, 1.0, None, False),
        wiring.Link(3, ("c", 1, 2), "A", "c", 4.10, 1.0, None, False),  # the same gene again
        wiring.Link(4, ("c", 9, 9), "D", "c", 4.05, 1.0, None, False),  # another element
    ]
    groups, distinct = wiring.keys_for(links, "element_kept", expression_classes=False)
    alts = wiring.alternatives(links, groups, distinct, 0.2)
    assert alts[0] == {"B"}
    assert alts[1] == {"A"}
    assert alts[2] == set()
    assert alts[3] == {"B"}
    assert alts[4] == set()


def test_balance_statistics():
    assert wiring.smd([0, 1, 2], [0, 1, 2]) == 0
    assert wiring.smd([0, 2], [1, 3]) == pytest.approx(1 / math.sqrt(2))
    assert wiring.smd([1], [1, 2]) is None
    assert wiring.ks([1, 2, 3], [1, 2, 3]) == 0
    assert wiring.ks([1, 2], [3, 4]) == 1
    assert wiring.ks([1, 2, 3, 4], [3, 4, 5, 6]) == pytest.approx(0.5)
    assert wiring.quantiles([0, 10], (0.5,)) == [5]


def test_expression_classes_keep_unknown_apart():
    links = [
        wiring.Link(0, ("c", 1, 2), "A", "c", 4.0, 1.0, None, True),
        wiring.Link(1, ("c", 1, 2), "B", "c", 4.0, 1.0, 2.0, False),
        wiring.Link(2, ("c", 1, 2), "C", "c", 4.0, 1.0, None, False),
    ]
    groups, _ = wiring.keys_for(links, "element_kept", expression_classes=True)
    assert groups[0] == groups[2] != groups[1]


def test_diagnostics_count_every_pair_and_hold_the_exact_checks():
    pairs = synthetic()
    expr = expression_for(pairs)
    a = wiring.assess(pairs, expr, "training", "element_kept", wiring.SCHEMES[2], draws=8)
    d = a["diagnostics"]
    m = d["meaningful"]
    assert m["pairs_evaluated"] == len(pairs)
    assert m["pairs_outside_population"] == sum(not wiring.in_population(p) for p in pairs)
    assert (
        m["immovable_links"]["all"] + round(m["movable_share"]["all"] * m["population_links"]["all"])
        == (m["population_links"]["all"])
    )
    assert all(d["exact"].values())
    assert d["balance"]["activity"]["regulated"]["smd"] in (0, None)  # the element is kept
    names = {c["check"] for c in a["verdict"]["checks"]}
    assert "regulated links moved per rewiring" in names
    assert "expression known share, regulated" in names
    # 30 regulated links cannot meet a floor of 100, so the synthetic set is not adequate
    assert not a["verdict"]["adequate"]
    assert any(c["check"] == "regulated links moved per rewiring" for c in a["verdict"]["failed"])


def test_adequacy_fails_a_missing_value_and_passes_only_when_all_met():
    pairs = synthetic(n_elements=60)
    a = wiring.assess(
        pairs, expression_for(pairs), "heldout_k562", "element_kept", wiring.SCHEMES[2], draws=4
    )
    diag = a["diagnostics"]
    diag["balance"]["distance"]["regulated"]["smd"] = None
    v = wiring.adequacy(diag, 1)
    assert not v["adequate"]
    assert any(c["check"] == "distance |SMD|, regulated" for c in v["failed"])
    loose = dict.fromkeys(wiring.THRESHOLDS, 0.0) | {
        "smd_max": math.inf,
        "ks_max": math.inf,
        "expression_unknown_gap_max": math.inf,
    }
    diag["balance"]["distance"]["regulated"]["smd"] = 0.0
    assert wiring.adequacy(diag, 0, loose)["adequate"]


def test_choose_takes_the_first_adequate_scheme_in_order():
    def fake(name, ok):
        return {"scheme": {"name": name}, "verdict": {"adequate": ok}}

    assert wiring.choose([fake("S1", False), fake("S2", True), fake("S3", True)]) == {"name": "S2"}
    assert wiring.choose([fake("S1", False), fake("S2", False)]) is None


def test_feasibility_never_scores_an_assignment(monkeypatch):
    def refuse(*a, **k):
        raise AssertionError("feasibility must not fit a model or compute a metric")

    for name in (
        "logistic_fit",
        "logistic_score",
        "benchmark_auprc",
        "average_precision",
        "auroc",
        "metrics",
    ):
        monkeypatch.setattr(crispri, name, refuse)
    pairs = synthetic()
    for null in wiring.NULLS:
        for scheme in wiring.SCHEMES:
            wiring.assess(pairs, expression_for(pairs), "training", null, scheme, draws=3)


def test_expression_by_symbol_leaves_unknown_genes_out():
    out = wiring.expression_by_symbol({"A": "ENSG1", "B": "ENSG2"}, {"ENSG1": math.e - 1})
    assert out == {"A": pytest.approx(1.0)}
