# SPDX-License-Identifier: AGPL-3.0-or-later
"""R8's coupling pretest: labels, shares, components, intervals and the training-only read."""

from __future__ import annotations

import gzip

import pytest

from genomeos.attribution import joint_pretest as jp
from genomeos.attribution import measured


def row(chrom, start, gene, outcome, drop):
    return jp.Row(chrom, start, start + 500, gene, "K562", outcome, drop)


def test_labels_follow_r2():
    assert jp.label_of(measured.DECREASE) is True
    assert jp.label_of(measured.NULL_INFORMATIVE) is False
    for o in (measured.INCREASE, measured.NULL_INCONCLUSIVE, measured.MISSING):
        assert jp.label_of(o) is None


def test_shares_read_partners_and_skip_missing():
    rows = [
        row("chr1", 100, "A", measured.DECREASE, 0.9),
        row("chr1", 100, "B", measured.INCREASE, 0.3),  # a partner whatever its outcome
        row("chr1", 100, "C", measured.DECREASE, None),  # missing: no share, not a partner
        row("chr1", 900, "A", measured.NULL_INFORMATIVE, 0.0),
    ]
    jp.couple(rows)
    a = rows[0].features
    assert a["element_competition"] == pytest.approx(0.9 / (1.2 + jp.TAU))
    assert a["gene_budget"] == pytest.approx(0.9 / (0.9 + jp.TAU))
    assert a["self_share"] == pytest.approx(0.9 / (0.9 + jp.TAU))
    assert a["element_partners"] == 1 and a["gene_partners"] == 1
    assert rows[2].features == {}


def test_components_join_through_shared_genes_and_perturbations():
    rows = [
        row("chr1", 100, "A", measured.DECREASE, 0.5),
        row("chr1", 900, "A", measured.DECREASE, 0.5),
        row("chr1", 900, "B", measured.DECREASE, 0.5),
        row("chr2", 100, "C", measured.DECREASE, 0.5),
    ]
    c = jp.components(rows)
    assert c[0] == c[1] == c[2] != c[3]


def test_gain_interval_sees_a_real_gain_and_not_a_null():
    labels = [i % 5 == 0 for i in range(400)]
    good = [1.0 if y else 0.0 for y in labels]
    noise = [((i * 7919) % 101) / 101 for i in range(400)]
    groups = [i // 10 for i in range(400)]
    g = jp.gain_interval(good, noise, labels, groups, n=200)
    assert g["interval"][0] > 0
    same = jp.gain_interval(noise, noise, labels, groups, n=200)
    assert same["gain"] == 0 and same["interval"] == [0, 0]


def test_evaluate_reads_no_coupling_where_none_was_planted():
    # labels follow the drop alone; the shares carry nothing beyond it
    rows = []
    for c in range(8):
        for e in range(30):
            for g in range(2):
                d = ((c * 31 + e * 7 + g * 3) % 17) / 10
                outcome = measured.DECREASE if d > 1.3 else measured.NULL_INFORMATIVE
                rows.append(
                    jp.Row(f"chr{c + 1}", e * 1000, e * 1000 + 500, f"G{c}_{e}_{g}", "K562", outcome, d)
                )
    res = jp.evaluate(rows, n=100)
    assert res["r8"] == "closes"
    assert res["reading"] in ("neither_passes", "void")


def test_evaluate_finds_planted_element_competition():
    # the positive is the element's stronger gene, and each element's scale varies tenfold, so the
    # absolute drop ranks poorly and the share within the perturbation ranks perfectly
    rows = []
    for c in range(10):
        for e in range(40):
            scale = 0.2 + ((c * 13 + e * 5) % 10) / 3
            for g, frac in enumerate((1.0, 0.4)):
                outcome = measured.DECREASE if g == 0 else measured.NULL_INFORMATIVE
                rows.append(
                    jp.Row(
                        f"chr{c + 1}",
                        e * 1000,
                        e * 1000 + 500,
                        f"G{c}_{e}_{g}",
                        "K562",
                        outcome,
                        scale * frac,
                    )
                )
    res = jp.evaluate(rows, n=200)
    assert res["verdict"]["element_competition"]["passes"]
    assert res["r8"] == "proceeds"


def test_score_refuses_a_heldout_pair():
    p = measured.CrispriPair(
        "chr1", 1, 2, "A", "K562", "d", "r", True, True, -0.5, 0.01, split=measured.HELDOUT
    )
    with pytest.raises(ValueError):
        jp.score([p])


def test_training_load_never_opens_the_heldout_file(monkeypatch):
    names = dict(measured.CRISPRI_SPLIT_OF)
    if not any((measured.CRISPRI_KNOWLEDGE / n).exists() for n in names):
        pytest.skip("the CRISPRi benchmark is not cached on this machine")
    opened = []
    real = gzip.open

    def spy(path, *a, **k):
        opened.append(str(path))
        return real(path, *a, **k)

    monkeypatch.setattr(measured.gzip, "open", spy)
    pairs = jp.load_training()
    assert pairs and all(p.split == measured.TRAINING for p in pairs)
    held = [n for n, s in names.items() if s == measured.HELDOUT]
    assert not any(h in o for h in held for o in opened)
