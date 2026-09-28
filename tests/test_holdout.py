# SPDX-License-Identifier: AGPL-3.0-or-later
"""The leave-one-source-out harness (item 13 C4): metrics, loci, the split, the leak guard, the
labellings and the call the pilot makes, on synthetic units only (no data/knowledge needed)."""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from genomeos.attribution import holdout as ho  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402


def _ap_brute(scores, y):
    """Average precision over distinct thresholds, by definition."""
    pairs = sorted(zip(scores, y, strict=True), key=lambda t: -t[0])
    n_pos = sum(y)
    ap, tp, fp, i = 0.0, 0, 0, 0
    while i < len(pairs):
        j = i
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            j += 1
        dtp = sum(t[1] for t in pairs[i:j])
        tp += dtp
        fp += (j - i) - dtp
        ap += (dtp / n_pos) * (tp / (tp + fp))
        i = j
    return ap


def _auc_brute(scores, y):
    pos = [s for s, t in zip(scores, y, strict=True) if t]
    neg = [s for s, t in zip(scores, y, strict=True) if not t]
    return sum((p > n) + 0.5 * (p == n) for p, n in itertools.product(pos, neg)) / (len(pos) * len(neg))


def test_average_precision_and_auroc_match_their_definitions_with_ties():
    scores = [0.9, 0.9, 0.5, 0.5, 0.5, 0.2, 0.1, 0.1]
    y = [1, 0, 1, 1, 0, 0, 1, 0]
    m = ho.binary_metrics(scores, y, np.ones((1, len(y))))
    assert m["average_precision"][0] == pytest.approx(_ap_brute(scores, y))
    assert m["auroc"][0] == pytest.approx(_auc_brute(scores, y))


def test_a_weight_of_two_is_the_unit_counted_twice():
    scores = [3.0, 2.0, 2.0, 1.0, 0.5]
    y = [1, 0, 1, 0, 1]
    w = np.array([[2.0, 1.0, 1.0, 3.0, 1.0]])
    rep = [s for s, k in zip(scores, w[0], strict=True) for _ in range(int(k))]
    ry = [t for t, k in zip(y, w[0], strict=True) for _ in range(int(k))]
    m = ho.binary_metrics(scores, y, w)
    assert m["average_precision"][0] == pytest.approx(_ap_brute(rep, ry))
    assert m["auroc"][0] == pytest.approx(_auc_brute(rep, ry))


def test_abstentions_rank_below_every_stated_score():
    # an abstention is never promoted above a stated zero
    m = ho.binary_metrics([None, 0.0, None], [1, 0, 0], np.ones((1, 3)))
    assert m["auroc"][0] == pytest.approx(0.25)


def test_spearman_is_the_rank_correlation():
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    assert ho.spearman(x, [2.0, 1.0, 4.0, 3.0, 5.0], np.ones((1, 5)))[0] == pytest.approx(0.8)
    assert ho.spearman(x, [5.0, 4.0, 3.0, 2.0, 1.0], np.ones((1, 5)))[0] == pytest.approx(-1.0)


def _u(source, chrom="chr1", start=0, end=100, **kw):
    return ho.Unit(source=source, chrom=chrom, start=start, end=end, **kw)


def test_loci_join_by_bin_and_for_pairs_by_gene():
    a = _u("crispri:X", start=100, end=200, gene="G1")
    b = _u("crispri:X", start=5_000_100, end=5_000_200, gene="G1")  # another bin, same gene
    c = _u("crispri:X", start=9_000_100, end=9_000_200, gene="G2")
    assert ho.loci([a, b, c]) == [0, 0, 1]
    d = _u("vista", start=150, end=250)
    e = _u("vista", start=5_000_100, end=5_000_200)
    assert ho.loci([d, e]) == [0, 1]  # element units never join on a gene


def _crispri(study, start, gene, outcome, split=ms.TRAINING, chrom="chr1", cell="K562", tss=None):
    return _u(
        f"crispri:{study}",
        chrom=chrom,
        start=start,
        end=start + 500,
        gene=gene,
        gene_id=f"ID{gene}",
        outcome=outcome,
        value=-0.3,
        cell=cell,
        split=split,
        study=study,
        tss=start + 10_000 if tss is None else tss,
    )


REFS = {
    "crispri:A": frozenset({"a 2019"}),
    "crispri:B": frozenset({"b 2020"}),
    "crispri:C": frozenset({"c 2021"}),
    "crispri:D": frozenset({"b 2020"}),
}


def _world():
    return {
        "crispri:A": (
            _crispri("A", 1000, "G1", ms.DECREASE),
            _crispri("A", 50_000, "G2", ms.NULL_INFORMATIVE),
        ),
        "crispri:B": (_crispri("B", 1200, "G1", ms.DECREASE), _crispri("B", 90_000, "G3", ms.DECREASE)),
        "crispri:C": (_crispri("C", 70_000, "G4", ms.DECREASE, split=ms.HELDOUT),),
        "crispri:D": (_crispri("D", 200_000, "G5", ms.DECREASE),),
        "lentimpra:K562": (_u("lentimpra:K562", start=1000, end=1200, value=2.0, cell="K562"),),
        "vista": (_u("vista", start=40_000, end=41_000, outcome="positive"),),
    }


def test_the_view_removes_the_source_its_siblings_and_masks_related_intervals_of_its_family():
    v = ho.evidence("crispri:A", _world(), REFS)
    assert "crispri:A" not in v.sources
    # B shares a base with A's first pair: masked; B's other pair stays
    assert [u.start for u in v.sources["crispri:B"]] == [90_000]
    assert v.masked == {"crispri:B": 1}
    # the held-out file is never evidence, whatever is held out
    assert "crispri:C" not in v.sources and v.heldout_file_pairs_dropped == 1
    # another family at the same interval is evidence in full
    assert len(v.sources["lentimpra:K562"]) == 1
    # a provenance sibling (same Reference) goes with the held-out study
    vb = ho.evidence("crispri:B", _world(), REFS)
    assert "crispri:D" not in vb.sources and "crispri:D" in vb.removed


def test_the_leak_guard_refuses_what_the_split_forbids():
    w = _world()
    ok = ho.rest_labels(ho.evidence("crispri:A", w, REFS))
    ho.check_provenance(ok, "crispri:A", REFS)
    with pytest.raises(ho.LeakError):
        ho.check_provenance(ho.Labels("x", lambda u, e: 0.0, frozenset({"crispri:A"})), "crispri:A", REFS)
    with pytest.raises(ho.LeakError):  # a same-family source read outside the masked view
        ho.check_provenance(ho.Labels("x", lambda u, e: 0.0, frozenset({"crispri:B"})), "crispri:A", REFS)
    with pytest.raises(ho.LeakError):
        ho.check_provenance(
            ho.Labels("x", lambda u, e: 0.0, frozenset({ho.CRISPRI_HELDOUT_FILE})), "vista", REFS
        )
    with pytest.raises(ho.LeakError):  # a view built for another held-out source
        ho.check_provenance(ho.rest_labels(ho.evidence("crispri:B", w, REFS)), "crispri:A", REFS)
    # another family is not a leak
    other = ho.Labels("x", lambda u, e: 0.0, frozenset({"lentimpra:K562", ho.MODEL}))
    ho.check_provenance(other, "crispri:A", REFS)


def test_endpoints_keep_the_five_outcomes_apart():
    dec = _crispri("A", 0, "G", ms.DECREASE)
    inc = _crispri("A", 0, "G", ms.INCREASE)
    null = _crispri("A", 0, "G", ms.NULL_INFORMATIVE)
    weak = _crispri("A", 0, "G", ms.NULL_INCONCLUSIVE)
    assert [ho.target(x, "decrease") for x in (dec, inc, null, weak)] == [1.0, None, 0.0, None]
    assert [ho.target(x, "increase") for x in (dec, inc, null, weak)] == [None, 1.0, 0.0, None]


def _big_source(n_loci=12, per=6, seed=1):
    rng = np.random.default_rng(seed)
    units = []
    for i in range(n_loci):
        for j in range(per):
            pos = i * 3_000_000 + j * 1000
            outcome = ms.DECREASE if rng.random() < 0.4 else ms.NULL_INFORMATIVE
            near = 1000 if outcome == ms.DECREASE else 90_000
            units.append(
                _crispri("S", pos, f"G{i}_{j}", outcome, tss=pos + near + int(rng.integers(0, 50_000)))
            )
    return {"crispri:S": tuple(units)}


def test_score_is_deterministic_cached_and_reports_coverage(tmp_path, monkeypatch):
    monkeypatch.setattr(ho, "CACHE_DIR", tmp_path)
    ho._MEMO.clear()
    world = _big_source()
    refs = {"crispri:S": frozenset({"s"})}
    lab = ho.distance_labels()
    a = ho.score(lab, "crispri:S", units=world, references=refs, resamples=200)
    ho._MEMO.clear()
    b = ho.score(lab, "crispri:S", units=world, references=refs, resamples=200)
    assert a == b and a["status"] == "scored" and a["loci"] == 12
    assert a["interval"][0] <= a["value"] <= a["interval"][1]
    assert a["coverage"]["share"] == 1.0 and a["benchmark_status"] == ho.STATUS
    assert any(Path(tmp_path, "score").iterdir())


def test_floors_describe_rather_than_score():
    world = _big_source(n_loci=5)
    r = ho.score(ho.distance_labels(), "crispri:S", units=world, references={}, resamples=50, cache=False)
    assert r["status"] == "described_not_scored" and ("loci" in r["why"] or "positives" in r["why"])


def test_compare_of_a_labelling_with_itself_is_zero():
    world = _big_source()
    lab = ho.distance_labels()
    other = ho.Labels("same", lab.predict, frozenset({"gencode_v50"}))
    c = ho.compare(lab, other, "crispri:S", units=world, references={}, resamples=100, cache=False)
    assert c["difference"] == 0 and c["interval"] == [0.0, 0.0]


def test_unchanged_labels_follow_direction_and_abstain_off_the_program():
    links = {
        "chr1": [
            ho.Link("chr1", 1000, 1300, "E1", "G1", "activates", 0.4, "K562"),
            ho.Link("chr1", 5000, 5300, "E2", "G2", "inhibits", 0.2, "HepG2"),
        ]
    }
    lab = ho.unchanged_labels(links)
    assert lab.model_dependent
    on = _crispri("A", 1100, "G1", ms.DECREASE)
    assert lab.predict(on, "decrease") == 0.4
    assert lab.predict(on, "increase") == 0.0  # the link is an activation
    assert lab.predict(_crispri("A", 1100, "OTHER", ms.DECREASE), "decrease") == 0.0
    assert lab.predict(_crispri("A", 20_000, "G1", ms.DECREASE), "decrease") is None
    assert lab.predict(_u("vista", start=5100, end=5200), "positive") == 0.2


def test_rest_orders_gene_support_then_element_support_then_distance(monkeypatch):
    monkeypatch.setattr(ho, "nearest_coding_distance", lambda c, p: 0)
    world = {
        "crispri:T": (_crispri("T", 1000, "G1", ms.DECREASE),),
        "gtex": (
            _u("gtex", start=1000, end=1400, gene="G1", gene_id="IDG1", outcome="associated", tss=11_000),
        ),
        "lentimpra:K562": (_u("lentimpra:K562", start=30_000, end=30_200, value=3.0, cell="K562"),),
        "crispri:H": (
            _crispri("H", 1100, "G1", ms.DECREASE),
            _crispri("H", 30_050, "G9", ms.DECREASE),
            _crispri("H", 60_000, "G8", ms.DECREASE),
        ),
    }
    refs = {k: frozenset({k}) for k in world if k.startswith("crispri")}
    view = ho.evidence("crispri:H", world, refs)
    assert "crispri:T" not in view.sources  # T shares bases with H's first pair: masked
    lab = ho.rest_labels(view)
    gene_supported, element_supported, nothing = world["crispri:H"]
    s1 = lab.predict(gene_supported, "decrease")  # GTEx supports the gene
    s2 = lab.predict(element_supported, "decrease")  # an active reporter tile
    s3 = lab.predict(nothing, "decrease")
    assert s1 >= 10 > s2 >= 1 > s3 > 0


PROGRAM = """module human.noncoding.chr1

element E1 {
  class: enhancer
  locus: chr1:1000-1300
  targets: G1
  activity: activates_target
}
rule E1 activates G1 { strength: 0.4; when: cell_type = K562; evidence: predicted "x" }
element E1_measured {
  class: enhancer
  locus: chr1:1000-1300
  activity: activates_target
}
rule E1_measured activates G1 { strength: 0.3; when: cell_type = K562; evidence: experimental "y" }
"""


def test_the_program_parser_reads_predicted_links_only(tmp_path, monkeypatch):
    d = tmp_path / "compiled"
    d.mkdir()
    (d / "noncoding_chr1.bio").write_text(PROGRAM)
    monkeypatch.setattr(ho, "CACHE_DIR", tmp_path / "cache")
    links = ho.compiled_links(d)
    assert [(x.element, x.gene, x.action, x.strength, x.cell) for x in links["chr1"]] == [
        ("E1", "G1", "activates", 0.4, "K562")
    ]
    kinds = [(k, n) for k, n, _ in ho.program_blocks(d / "noncoding_chr1.bio")]
    assert kinds == [("element", "E1"), ("rule", "E1"), ("element", "E1_measured"), ("rule", "E1_measured")]


def test_the_registration_is_complete_and_says_what_it_may_not_be_called():
    r = ho.registration()
    assert r["status"] == "internal development benchmark"
    assert any("external" in x for x in r["may_not_be_called"])
    assert set(r["labellings"]) == {"rest", "distance", "unchanged"}
    assert r["resamples"] == 1000 and r["level"] == 0.95 and not math.isnan(r["seed"])
