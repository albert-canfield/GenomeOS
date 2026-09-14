"""The epigenome layer, offline: file choice, the bigBed reader, bins, records, states, stats."""

from __future__ import annotations

import gzip
import math
import struct
import zlib
from array import array

import pytest

from genomeos.genome import epigenome as ep
from genomeos.genome import reader


def _f(acc, output, fmt="bed", ftype="bed narrowPeak", reps=(1, 2), size=100, default=False):
    return {
        "accession": acc,
        "output_type": output,
        "file_format": fmt,
        "file_type": ftype,
        "file_size": size,
        "biological_replicates": list(reps),
        "preferred_default": default,
        "cloud_metadata": {"url": f"https://example.org/{acc}"},
    }


def test_histone_files_default_analysis_replicated_and_pooled():
    files = [
        _f("P1", "pseudoreplicated peaks", size=900, default=True),
        _f("P2", "replicated peaks", size=100),
        _f("P3", "replicated peaks", size=500),
        _f("OLD", "replicated peaks", size=9999),
        _f("S1", "fold change over control", fmt="bigWig", ftype="bigWig", reps=(1,), size=900),
        _f("S2", "fold change over control", fmt="bigWig", ftype="bigWig", reps=(1, 2), size=100),
    ]
    got = ep.choose_histone_files(files, {"P1", "P2", "P3", "S1", "S2"})
    assert got["peaks"]["accession"] == "P3"  # replicated, larger; OLD is outside the analysis
    assert got["signal"]["accession"] == "S2"  # pooled replicates beat size
    assert got["signal"]["url"].endswith("S2")


def test_treatment_rule_keeps_differentiation():
    assert ep._treated("Homo sapiens SK-N-SH treated with 6 uM all-trans-retinoic acid")
    assert not ep._treated("Homo sapiens cardiac muscle cell originated from H7 treated with BMP4")
    assert not ep._treated("Homo sapiens K562")


def _opt(mark, biosample, acc, peaks=True, signal=True, size=10):
    return {
        "mark": mark,
        "experiment": acc,
        "biosample": biosample,
        "peaks": {"accession": acc + "p", "size": size} if peaks else None,
        "signal": {"accession": acc + "s"} if signal else None,
    }


def test_pick_marks_prefers_one_biosample():
    rues, h7 = "cardiac muscle cell originated from RUES2", "cardiac muscle cell originated from H7"
    options = {
        "H3K4me3": [_opt("H3K4me3", h7, "A", size=99), _opt("H3K4me3", rues, "B")],
        "H3K27ac": [_opt("H3K27ac", rues, "C")],
        "H3K27me3": [_opt("H3K27me3", h7, "D", size=99), _opt("H3K27me3", rues, "E")],
        "H3K9me3": [],
        "H3K4me1": [_opt("H3K4me1", h7, "F", peaks=False, signal=False)],
    }
    got = ep.pick_marks(options)
    assert got["H3K4me3"]["experiment"] == "B"
    assert got["H3K27me3"]["experiment"] == "E"
    assert got["H3K9me3"]["value"] == ep.UNKNOWN
    assert got["H3K4me1"]["value"] == ep.UNKNOWN and "no GRCh38 peaks" in got["H3K4me1"]["reason"]


def _bigbed(path, chrom, length, rows):
    """A minimal little-endian bigBed: one chromosome, one leaf, one zlib data block."""
    block = b"".join(struct.pack("<III", 0, s, e) + rest.encode() + b"\0" for s, e, rest in rows)
    comp = zlib.compress(block)
    key_size = 8
    chrom_tree = struct.pack("<IIIIQQ", 0x78CA8C91, 1, key_size, 8, 1, 0)
    chrom_tree += (
        struct.pack("<BBH", 1, 0, 1) + chrom.encode().ljust(key_size, b"\0") + struct.pack("<II", 0, length)
    )
    data_off = 64 + len(chrom_tree)
    data = struct.pack("<Q", len(rows)) + comp
    index_off = data_off + len(data)
    rtree = struct.pack("<IIQIIIIQII", 0x2468ACE0, 1, 1, 0, rows[0][0], 0, rows[-1][1], index_off, 1, 0)
    rtree += struct.pack("<BBH", 1, 0, 1) + struct.pack(
        "<IIIIQQ", 0, rows[0][0], 0, rows[-1][1], data_off + 8, len(comp)
    )
    header = struct.pack(
        "<IHHQQQHHQQIQ", ep.BIGBED_MAGIC, 4, 0, 64, data_off, index_off, 14, 9, 0, 0, len(block) + 1, 0
    )
    path.write_bytes(header + chrom_tree + data + rtree)


def test_bigbed_reader_and_bedmethyl_bins(tmp_path):
    rows = [
        (100, 101, '"0"\t10\t+\t100\t101\t0,255,0\t10\t80\tCG\tCG\t40'),
        (101, 102, '"0"\t4\t-\t101\t102\t0,255,0\t4\t0\tCG\tCG\t40'),  # below MIN_COVERAGE
        (450, 451, '"0"\t20\t+\t450\t451\t0,255,0\t20\t50\tCG\tCG\t40'),
    ]
    p = tmp_path / "t.bb"
    _bigbed(p, "chr21", 1000, rows)
    bb = ep.BigBed(p)
    try:
        got = list(bb.rows("chr21", [(0, 1000)]))
        assert list(bb.rows("chr9", [(0, 1000)])) == []
    finally:
        bb.close()
    assert [(s, e) for s, e, _ in got] == [(100, 101), (101, 102), (450, 451)]
    calls = [(s, *ep.bedmethyl_call(r)) for s, _e, r in got]
    assert calls[0] == (100, 10, 0.8)
    cols = ep.bin_calls(calls, 1000)
    m = ep.methylation_over(cols, 0, 200)
    assert m == {"fraction": 0.8, "cpg_calls": 2, "cpg_calls_covered": 1, "mean_coverage": 10.0}
    both = ep.methylation_over(cols, 0, 1000)
    assert both["fraction"] == 0.65 and both["cpg_calls_covered"] == 2
    assert ep.methylation_over(cols, 600, 800)["fraction"] is None
    bad = tmp_path / "bad.bb"
    bad.write_bytes(b"\0" * 64)
    with pytest.raises(ValueError):
        ep.BigBed(bad)


def test_profile_mean_skips_bins_without_value():
    prof = array("f", [1.0, math.nan, 3.0, 5.0])
    assert ep.profile_mean(prof, 0, 400) == 1.0
    assert ep.profile_mean(prof, 0, 600) == 2.0
    assert ep.profile_mean(prof, 250, 260) is None
    assert ep.profile_mean(prof, 450, 460) == 3.0
    assert ep.profile_mean(array("f", [math.nan]), 0, 100) is None


@pytest.mark.parametrize(
    ("called", "state"),
    [
        ({"H3K4me3": True, "H3K27me3": True}, "bivalent"),
        ({"H3K4me3": True, "H3K27ac": True}, "active promoter"),
        ({"H3K27ac": True, "H3K4me1": True}, "active enhancer"),
        ({"H3K4me1": True, "H3K27me3": True}, "poised enhancer"),
        ({"H3K4me1": True}, "primed enhancer"),
        ({"H3K27me3": True}, "Polycomb-repressed"),
        ({"H3K9me3": True}, "heterochromatin"),
        ({}, "no mark"),
        ({"H3K9me3": None}, ep.UNKNOWN),
    ],
)
def test_chromatin_state_rules(called, state):
    full = {m: False for m in ep.MARKS} | called
    assert ep.chromatin_state(full) == state


def _write_peaks(path, rows, header="# test\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write(header)
        for s, e, v in rows:
            fh.write(f"{s}\t{e}\t{v}\n")


def test_record_fields_evidence_and_unknown(tmp_path, monkeypatch):
    monkeypatch.setattr(ep, "CACHE", tmp_path / "knowledge")
    monkeypatch.setattr(ep, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(reader, "RESULTS", tmp_path / "results")
    _write_peaks(tmp_path / "results" / "dnase_K562_chr21.bed.gz", [(900, 1100, 40.0)], "# ENCFFDNASE K562\n")
    for mark in ("H3K4me3", "H3K27ac", "H3K4me1", "H3K9me3"):
        _write_peaks(ep.peaks_path("K562", mark, "chr21"), [(950, 1050, 7.0)] if mark != "H3K9me3" else [])
    prof = array("f", [0.5] * 10)
    sp = ep.signal_path("K562", "H3K27ac", "chr21")
    sp.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(sp, "wb") as fh:
        fh.write(prof.tobytes())
    entry = {
        "experiment": "ENCSRX",
        "biosample": "Homo sapiens K562",
        "peaks": {"accession": "ENCFFP", "output_type": "replicated peaks"},
        "signal": {"accession": "ENCFFS", "replicates": [1, 2]},
    }
    manifest = {
        "cell_types": {
            "K562": {
                "marks": {m: dict(entry) for m in ep.MARKS},
                "methylation": {"value": ep.UNKNOWN, "reason": "no released GRCh38 WGBS"},
            }
        }
    }
    layer = ep.Layer("chr21", ["K562", "keratinocyte"], manifest=manifest)
    rec_k, rec_kera = ep.epigenome_at("chr21", 1000, 1200, ["K562", "keratinocyte"], layer=layer)
    assert rec_k["openness"]["value"] == "open" and rec_k["openness"]["source"] == "ENCFFDNASE"
    k27ac = rec_k["marks"]["H3K27ac"]
    assert k27ac["peak"] is True and k27ac["fold_change"] == 0.5 and k27ac["confidence"] == 0.9
    assert k27ac["evidence"].startswith("experimental")
    assert rec_k["marks"]["H3K9me3"]["peak"] is False
    # H3K27me3 has an experiment but this chromosome was never read: UNKNOWN, not absent
    assert rec_k["marks"]["H3K27me3"]["value"] == ep.UNKNOWN
    assert rec_k["state"]["value"] == "active promoter"
    assert rec_k["methylation"] == {"value": ep.UNKNOWN, "reason": "no released GRCh38 WGBS"}
    assert rec_kera["openness"]["value"] == ep.UNKNOWN
    assert all(v["value"] == ep.UNKNOWN for v in rec_kera["marks"].values())

    # the chromosome summary over the same cache: one read promoter, one silent, two elements
    from types import SimpleNamespace

    from genomeos.coords import Locus, Strand
    from genomeos.genome.regulatory import CCRE

    genes = {
        "A": SimpleNamespace(
            symbol="A", type="protein_coding", locus=Locus("chr21", 1000, 5000, Strand.PLUS)
        ),
        "B": SimpleNamespace(
            symbol="B", type="protein_coding", locus=Locus("chr21", 7000, 9000, Strand.MINUS)
        ),
    }
    ccres = [CCRE("chr21", 960, 1040, "E1", "PLS", False), CCRE("chr21", 5000, 5100, "E2", "dELS", True)]
    s = ep.summarise_chromosome(
        "chr21", SimpleNamespace(genes=genes), ccres, cell_types=["K562"], manifest=manifest
    )
    row = s["cell_types"]["K562"]
    assert row["promoters"]["read"]["n"] == 1 and row["promoters"]["silent"]["n"] == 1
    assert row["promoters"]["read"]["H3K4me3"] == {"n": 1, "share": 1.0}
    assert row["promoters"]["silent"]["H3K4me3"] == {"n": 0, "share": 0.0}  # measured zero, not missing
    assert set(row["promoters"]) == {"read", "silent"}  # no empty "unread" bucket
    assert "H3K27me3" in row["marks_unknown"] and row["methylation"] is False
    assert row["registry_classes"]["PLS"]["open_share"] == 1.0


def test_coverage_table_names_what_is_missing():
    cells = {
        "A": {"dnase": "X", "marks": {m: {"experiment": "E"} for m in ep.MARKS}, "methylation": {"cpg": {}}},
        "B": {
            "dnase": "Y",
            "marks": {m: {"experiment": "E"} for m in ep.MARKS},
            "methylation": {"value": ep.UNKNOWN},
        },
    }
    t = ep.coverage_table(cells)
    assert t["complete"] == ["A"] and t["missing"] == {"B": ["methylation"]}
    assert t["with_all_five_marks"] == 2 and t["with_methylation"] == 1


def test_small_statistics():
    assert ep.auc([0.1, 0.2, 0.8, 0.9], [0, 0, 1, 1]) == 1.0
    assert ep.auc([0.5, 0.5], [0, 1]) == 0.5
    assert ep.auc([1, 2], [1, 1]) is None
    assert ep.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    x = [[1.0, float(i)] for i in range(20)]
    y = [2.0 + 3.0 * i for i in range(20)]
    w = ep.ridge_fit(x, y, lam=0.0)
    assert w == pytest.approx([2.0, 3.0], abs=1e-6)
    pred = ep.cross_validated(x, y, list(range(20)), folds=4, lam=0.0)
    assert pred == pytest.approx(y, abs=1e-6)
    z = ep.standardise([[1.0], [3.0]])
    assert z == [[1.0, -1.0], [1.0, 1.0]]
    gc, cpg = ep.gc_cpg("ACGTNNCG")
    assert gc == pytest.approx(4 / 6) and cpg == pytest.approx(200 / 6)
    assert ep.stratum(0.42, 1.5) == "gc8_cpg2"
