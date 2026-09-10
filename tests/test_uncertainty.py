from genomeos.ir import Evidence, EvidenceKind, Parameter
from genomeos.lang import parse_file
from genomeos.runtime.uncertainty import UncertaintyReport, report_for_ageing, report_for_network


def test_levels_and_unknown():
    m = parse_file("data/demo/repressilator.bio")
    rep = report_for_network(m, m.rules)
    d = rep.to_dict()
    assert d["molecular"]["label"] in ("high", "medium") and d["molecular"]["items"] == 9
    assert d["cellular"]["label"] == "UNKNOWN" and d["organism"]["label"] == "UNKNOWN"
    assert "UNKNOWN" in rep.format()


def test_weights_and_weakest():
    rep = UncertaintyReport()
    rep.add("cellular", 0.9, EvidenceKind.EXPERIMENTAL, "good")
    rep.add("cellular", 0.9, EvidenceKind.NONE, "guess")
    assert rep.levels["cellular"].weakest == "guess"
    assert abs(rep.levels["cellular"].confidence - 0.45) < 1e-9


def test_ageing_report_discounts_tissue():
    params = [Parameter("a", 1, "", Evidence(EvidenceKind.EXPERIMENTAL, "x"), 0.8)]
    rep = report_for_ageing(params)
    assert rep.levels["cellular"].confidence > rep.levels["tissue"].confidence
    assert rep.levels["organism"].confidence is None
