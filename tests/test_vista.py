# SPDX-License-Identifier: AGPL-3.0-or-later
"""VISTA loci parsed, tissues grouped, the registry and node read blind, positives against negatives."""

import gzip
import io
from types import SimpleNamespace

from genomeos.attribution import vista
from genomeos.coords import Locus, Strand

LOCI = (
    '"vista_id"\t"curation_status"\t"tissue"\t"experiments"\t"backbone"\t"stage"\t"assembly"\t"coord"'
    '\t"coordinate_hg38"\t"coordinate_mm10"\n'
    '"hs1"\t"positive"\t"cn;hb;lb;nt"\t1\t"hZR"\t"e11.5"\t"hg38"\t"chr16:86396481-86397120"'
    '\t"chr16:86396481-86397120"\t"chr8:121002895-121003512"\n'
    '"hs2"\t"negative"\t""\t1\t"hZR"\t"e11.5"\t"hg38"\t"chr16:85586489-85588130"'
    '\t"chr16:85586489-85588130"\t"chr8:120516455-120517853"\n'
    '"mm3"\t"positive"\t"ht"\t1\t"hZR"\t"e11.5"\t"mm10"\t"chr1:1-100"\t""\t"chr1:1-100"\n'
)


def test_parse_loci_keeps_human_hg38_only():
    els = vista.parse_loci(io.StringIO(LOCI))
    assert [e.id for e in els] == ["hs2", "hs1"]  # sorted by position
    hs1 = els[1]
    assert (hs1.chrom, hs1.start, hs1.end) == ("chr16", 86396480, 86397120)
    assert hs1.status == "positive" and hs1.tissues == ("cn", "hb", "lb", "nt")
    assert hs1.groups == {"neural", "limb and mesenchyme"}
    assert els[0].tissues == () and els[0].groups == set()
    assert vista.parse_loci(io.StringIO(LOCI), "chr1") == []


def test_load_loci_reads_the_cached_gzip(tmp_path):
    (tmp_path / "locus.tsv.gz").write_bytes(gzip.compress(LOCI.encode()))
    assert [e.id for e in vista.load_loci("chr16", tmp_path)] == ["hs2", "hs1"]


def test_track_groups():
    assert vista.group_of_track("Brain_Cerebellum") == "neural"
    assert vista.group_of_track("Heart_Left_Ventricle") == "heart"
    assert vista.group_of_track("chondrocyte") == "limb and mesenchyme"
    assert vista.group_of_track("K562") is None and vista.group_of_track(None) is None


def _gene(symbol, start, end, strand=Strand.PLUS, kind="protein_coding"):
    return SimpleNamespace(symbol=symbol, type=kind, locus=Locus("chr9", start, end, strand))


def _ctx():
    ccres = [
        SimpleNamespace(chrom="chr9", start=1_000, end=1_300, id="E1", cls="dELS", ctcf_bound=False),
        SimpleNamespace(chrom="chr9", start=5_000, end=5_200, id="E2", cls="CTCF-only", ctcf_bound=True),
    ]
    domains = [SimpleNamespace(id="chr9:D1", start=0, end=10_000, genes=["A", "B"])]
    genes = {
        "A": _gene("A", 3_000, 4_000),
        "B": _gene("B", 8_000, 9_000, Strand.MINUS),
        "L": _gene("L", 1, 2, kind="lncRNA"),
    }
    el = SimpleNamespace(id="E1", targets=[{"gene": "A", "distance": 1850, "basis": "nearest TSS in domain"}])
    return SimpleNamespace(
        chrom="chr9",
        ccres=ccres,
        domains=domains,
        annotation=SimpleNamespace(genes=genes),
        element_by_id={"E1": el},
    )


def test_annotate_reads_registry_node_and_constraint_blind():
    ctx = _ctx()
    els = [
        vista.VistaElement("hs10", "chr9", 900, 1_500, "positive", ("fb",), 1),
        vista.VistaElement("hs11", "chr9", 5_100, 5_600, "negative", (), 1),
        vista.VistaElement("hs12", "chr9", 20_000, 20_500, "negative", (), 1),
    ]
    stats = [
        SimpleNamespace(bases=600, fraction_above=0.4, mean=1.2),
        SimpleNamespace(bases=500, fraction_above=0.05, mean=0.1),
        SimpleNamespace(bases=0, fraction_above=None, mean=0.0),
    ]
    rows = vista.annotate(ctx, els, stats)
    r0, r1, r2 = rows
    assert r0["enhancer_like"] and r0["ccre_classes"] == ["dELS"] and r0["inferred"]["gene"] == "A"
    assert r0["domain"] == "chr9:D1" and r0["constrained_fraction"] == 0.4 and r0["groups"] == ["neural"]
    assert not r1["enhancer_like"] and r1["any_ccre"] and r1["ccre_classes"] == ["CTCF-only"]
    assert (
        r1["inferred"]["gene"] == "A" and r1["inferred"]["distance"] == 2_350
    )  # nearest coding TSS in the node
    assert r2["domain"] == "" and r2["inferred"] is None and r2["constrained_fraction"] is None


def test_summarise_separates_and_gives_a_chance_level():
    rows = [
        {
            "status": "positive",
            "enhancer_like": True,
            "any_ccre": True,
            "domain": "d",
            "groups": ["neural"],
            "constrained_fraction": 0.5,
            "predicted": {"strength": "strong", "tissue": "brain"},
            "verdict_coding": "agrees with nearest TSS in domain",
            "predicted_group": "neural",
            "tissue_agrees": True,
        },
        {
            "status": "positive",
            "enhancer_like": True,
            "any_ccre": True,
            "domain": "d",
            "groups": ["heart"],
            "constrained_fraction": 0.3,
            "predicted": {"strength": "weak", "tissue": "heart"},
            "verdict_coding": "gene outside the domain",
            "predicted_group": "heart",
            "tissue_agrees": True,
        },
        {
            "status": "positive",
            "enhancer_like": False,
            "any_ccre": False,
            "domain": "",
            "groups": ["neural"],
            "constrained_fraction": 0.1,
            "predicted": None,
            "predicted_group": None,
            "tissue_agrees": None,
        },
        {
            "status": "negative",
            "enhancer_like": False,
            "any_ccre": True,
            "domain": "d",
            "groups": [],
            "constrained_fraction": 0.1,
            "predicted": None,
            "predicted_group": None,
            "tissue_agrees": None,
        },
        {
            "status": "negative",
            "enhancer_like": True,
            "any_ccre": True,
            "domain": "d",
            "groups": [],
            "constrained_fraction": None,
            "predicted": {"strength": "weak", "tissue": "K562"},
            "verdict_coding": "no predicted effect",
            "predicted_group": None,
            "tissue_agrees": None,
        },
    ]
    s = vista.summarise(rows)
    p, n = s["positive"], s["negative"]
    assert p["elements"] == 3 and n["elements"] == 2
    assert p["enhancer_like"] == 0.667 and n["enhancer_like"] == 0.5
    assert p["constrained"] == 0.667 and n["constrained"] == 0.0
    assert p["fraction_with_target"] == 0.667 and n["fraction_with_target"] == 0.5
    assert p["strong"] == 1 and p["coding_target_agrees_with_nearest"] == 0.5
    assert p["tissue_judged"] == 2 and p["tissue_agrees"] == 1.0
    # shuffled labels: the neural prediction matches the other neural positive (1 of 1 other-labelled set
    # differing from its own is {"heart"}: 0), the heart prediction matches none of the neural sets: 0
    assert p["tissue_agrees_expected"] == 0.0
    assert s["separation"]["enhancer_like"] == 0.167 and s["separation"]["constrained"] == 0.667
    assert s["positive_groups"] == {"neural": 2, "heart": 1}
