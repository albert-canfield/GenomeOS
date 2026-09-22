# SPDX-License-Identifier: AGPL-3.0-or-later
"""Two people's variants inside regulatory elements compared, with the element's target and prediction."""

from types import SimpleNamespace

from genomeos.attribution.eqtl import Intervals
from genomeos.genome import regdiff


def _elements(monkeypatch):
    e1 = SimpleNamespace(
        id="E1",
        cls="enhancer",
        locus=SimpleNamespace(start=1000, end=1300),
        targets=[{"gene": "APP", "basis": "nearest TSS in domain", "distance": 12000}],
    )
    e2 = SimpleNamespace(
        id="E2",
        cls="promoter",
        locus=SimpleNamespace(start=5000, end=5400),
        targets=[{"gene": "SOD1", "basis": "promoter of", "distance": 100}],
    )
    e3 = SimpleNamespace(id="E3", cls="insulator", locus=SimpleNamespace(start=9000, end=9200), targets=[])
    by_id = {"E1": e1, "E2": e2, "E3": e3}
    iv = Intervals()
    for e in by_id.values():
        iv.add("chr21", e.locus.start, e.locus.end, e.id)
    monkeypatch.setattr(regdiff, "chromosome_elements", lambda chrom, reference=None: (by_id, iv.freeze()))
    monkeypatch.setattr(regdiff, "_reference_base_reader", lambda chrom, reference=None: (None, None))


def test_regulatory_diff_places_variants_and_ranks(tmp_path, monkeypatch):
    _elements(monkeypatch)
    calls = {
        "A": {(1100, "C", "T"): "het", (5100, "G", "A"): "hom", (7000, "A", "G"): "het"},  # 7000: no element
        "B": {(1100, "C", "T"): "het", (9100, "T", "C"): "het"},
    }
    monkeypatch.setattr(regdiff, "vcf_path", lambda name, chrom, root: f"{name}.vcf")
    monkeypatch.setattr(regdiff, "_genotypes", lambda path, base_at=None: calls[path.split(".")[0]])
    import genomeos.predict.enhancer_target as et

    def fake_pred(chrom, eid, cache=None):
        return {"gene": "SOD1", "log2_fold_change": -0.4, "tissue": "liver"} if eid == "E2" else None

    monkeypatch.setattr(et, "cached_prediction", fake_pred)
    d = regdiff.regulatory_diff("A", "B", "chr21", tmp_path)
    assert d["variants_in_elements_a"] == 2 and d["variants_in_elements_b"] == 2 and d["shared"] == 1
    assert d["only_a"] == 1 and d["only_b"] == 1
    assert d["with_predicted_target_only_a"] == 1 and d["by_class"]["only_b"] == {"insulator": 1}
    assert d["genes"][0]["gene"] == "SOD1"  # predicted target and promoter rank first
    md = regdiff.render(d)
    assert (
        "+ 5,100 G>A (homozygous) in promoter E2 → SOD1 (promoter of, 100 bp)  "
        "[deletion moves SOD1 -0.40 in liver]" in md
    )
    assert "- 9,100 T>C (heterozygous) in insulator E3" in md and "no coding target in the node" in md


def test_add_constraint_marks_constrained_bases(monkeypatch):
    import genomeos.attribution.constraint as c

    stats = [SimpleNamespace(bases=1, mean=3.1), SimpleNamespace(bases=0, mean=0.0)]
    monkeypatch.setattr(
        c, "phylop_over_blocks", lambda chrom, intervals, threshold=None, progress=None: (stats, {"mb": 0})
    )
    rows = [{"pos": 10}, {"pos": 20}, {"pos": 10}]
    regdiff.add_constraint("chr21", rows)
    assert rows[0]["phylop"] == 3.1 and rows[2]["phylop"] == 3.1 and rows[1]["phylop"] is None
    assert regdiff._constrained(rows[0]) and not regdiff._constrained(rows[1])
