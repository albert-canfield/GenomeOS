"""Peptide/HLA binding predictors: optional, listed, never invented.

These tests must pass whether or not MHCflurry or NetMHCpan is installed, so
they exercise the contract (availability, the honest default, the parser, the
Answer shape) rather than the models.
"""

from __future__ import annotations

import pytest

from genomeos.therapeutics.binding import (
    STRONG_RANK,
    WEAK_RANK,
    Binding,
    MhcflurryPredictor,
    NetMHCpanPredictor,
    choose,
    parse_netmhcpan,
    status,
)
from genomeos.therapeutics.providers import NoNeoantigenPredictor

NETMHCPAN_OUTPUT = """\
# NetMHCpan version 4.1b

HLA-A*02:01 : Distance to training data  0.000

-----------------------------------------------------------------------------------
 Pos         MHC        Peptide  Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL BindLevel
-----------------------------------------------------------------------------------
   1 HLA-A*02:01      KLVFFAEDV KLVFFAEDV  0  0  0  0  0    KLVFFAEDV         PEPLIST 0.8462100    0.089 <= SB
   1 HLA-A*02:01      KLVFFAEDK KLVFFAEDK  0  0  0  0  0    KLVFFAEDK         PEPLIST 0.0123000    1.500 <= WB
   1 HLA-A*02:01      DDDDDDDDD DDDDDDDDD  0  0  0  0  0    DDDDDDDDD         PEPLIST 0.0000100   45.000
-----------------------------------------------------------------------------------
"""


def test_every_predictor_is_listed_with_its_licence_and_how_to_enable_it():
    st = status()
    names = [p["name"] for p in st["predictors"]]
    assert names == ["MHCflurry 2", "NetMHCpan 4.1"]
    for p in st["predictors"]:
        assert p["licence"] and p["how"] and "what" in p
        assert isinstance(p["available"], bool)
    # Apache-licensed one is installable from here; the academic one never is
    assert "Apache" in st["predictors"][0]["licence"]
    assert "academic" in st["predictors"][1]["licence"]
    assert "not ship" in st["predictors"][1]["licence"]
    assert "presented peptide" in st["caveat"]


def test_choosing_none_gives_the_honest_default():
    assert isinstance(choose("none"), NoNeoantigenPredictor)


def test_an_absent_predictor_answers_unavailable_with_the_reason_never_an_affinity():
    p = NetMHCpanPredictor(binary=None)
    assert p.available is False
    a = p.predict(["KLVFFAEDV"], ["HLA-A*02:01"])
    assert a.available is False and a.data is None
    assert "not installed" in a.reason and "healthtech.dtu.dk" in a.reason


def test_a_predictor_with_no_peptides_or_no_alleles_says_so(monkeypatch):
    p = NetMHCpanPredictor(binary=__file__)  # any existing path counts as installed
    assert p.available is True
    assert "no peptides" in p.predict([], ["HLA-A*02:01"]).reason
    assert "no peptides" in p.predict(["KLVFFAEDV"], []).reason


def test_netmhcpan_output_is_parsed_by_header_not_by_column_position():
    b = parse_netmhcpan(NETMHCPAN_OUTPUT)
    assert len(b) == 3
    strong, weak, junk = b
    assert strong.peptide == "KLVFFAEDV" and strong.allele == "HLA-A*02:01"
    assert strong.percentile_rank == pytest.approx(0.089)
    assert strong.binder and strong.strong
    assert weak.binder and not weak.strong  # 1.5 is a binder, not a strong one
    assert not junk.binder  # 45.0 is not a binder


def test_thresholds_are_the_documented_class_one_conventions():
    assert WEAK_RANK == 2.0 and STRONG_RANK == 0.5
    assert Binding("P", "HLA-A*02:01", 2.0).binder
    assert not Binding("P", "HLA-A*02:01", 2.01).binder


def test_mhcflurry_reports_the_model_download_when_it_cannot_run(monkeypatch):
    p = MhcflurryPredictor()
    monkeypatch.setattr(MhcflurryPredictor, "installed", staticmethod(lambda: True))
    monkeypatch.setattr(
        MhcflurryPredictor, "_load", lambda self: (_ for _ in ()).throw(RuntimeError("Missing file"))
    )
    a = p.predict(["KLVFFAEDV"], ["HLA-A*02:01"])
    assert a.available is False and "mhcflurry-downloads fetch" in a.reason


def test_a_working_predictor_returns_binders_and_computational_evidence(monkeypatch):
    class FakeFrame:
        @staticmethod
        def to_dict(_orient):
            return [
                {"peptide": "KLVFFAEDV", "best_allele": "HLA-A*02:01", "presentation_percentile": 0.2},
                {"peptide": "DDDDDDDDD", "best_allele": "HLA-A*02:01", "presentation_percentile": 60.0},
            ]

    class FakeModel:
        @staticmethod
        def predict(peptides, alleles, verbose=0):
            return FakeFrame()

    p = MhcflurryPredictor()
    monkeypatch.setattr(MhcflurryPredictor, "installed", staticmethod(lambda: True))
    monkeypatch.setattr(MhcflurryPredictor, "_load", lambda self: FakeModel())
    a = p.predict(["KLVFFAEDV", "DDDDDDDDD"], ["HLA-A*02:01"])
    assert a.available is True
    assert a.data["predictor"] == "MHCflurry 2"
    assert len(a.data["binders"]) == 1 and a.data["strong_binders"] == 1
    assert a.data["best"]["peptide"] == "KLVFFAEDV"
    assert a.evidence and a.evidence[0].level == "computational"
    assert a.evidence[0].source_type == "prediction"
    assert "presented peptide" in a.data["caveat"]
    assert a.data["processing_modelled"] is True  # a presentation percentile includes processing


def test_a_binding_only_score_does_not_claim_antigen_processing(monkeypatch):
    """NetMHCpan's rank is binding; only a presentation score may claim processing."""
    from genomeos.therapeutics.binding import _summary

    d = _summary("NetMHCpan 4.1", [Binding("KLVFFAEDV", "HLA-A*02:01", 0.1)], ["HLA-A*02:01"])
    assert d["processing_modelled"] is False
    assert "not modelled here" in d["caveat"]
