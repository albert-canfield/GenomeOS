# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, stage 3: the constraint-dose view on a second source.

Module: genomeos/attribution/simpler_dose.py.

Synthetic checks of each piece the registration names: the stage call, the coverage-TPM and its strand
check, the endpoint, the matched-distance control on stages, the dose arms and their census, and the
whole score on planted data, including the outcomes that withhold a reading. The element reader is run
on tables that lack the source's own element-to-gene columns, so it cannot be reading them. The local
source is checked only when the git-ignored copies are present.
"""

from __future__ import annotations

import json
import struct
import zlib
from pathlib import Path

import numpy as np
import pytest

from genomeos.attribution import simpler
from genomeos.attribution import simpler_dose as sd
from genomeos.attribution import simpler_tiers as st
from genomeos.attribution.bigwig import LeafItem

ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------------------------- elements


def test_stage_of_takes_the_highest_and_ties_go_earlier() -> None:
    assert sd.stage_of([1.0, 5.0, 2.0, 0.0, 0.0, 4.9]) == "l1"
    assert sd.stage_of([3.0, 1.0, 3.0, 0.0, 0.0, 0.0]) == "emb"
    assert sd.stage_of([0.0] * 5 + [0.1]) == "ya"


def _write(path: Path, head: list[str], rows: list[list[str]]) -> None:
    path.write_text("\t".join(head) + "\n" + "".join("\t".join(r) + "\n" for r in rows))


def test_read_elements_needs_no_gene_column_and_keeps_dynamic_only(tmp_path: Path) -> None:
    fig1, fig2, fig4 = tmp_path / "f1.txt", tmp_path / "f2.txt", tmp_path / "f4.txt"
    heights = [f"atac_wt_{s}_height" for s in sd.STAGES]
    _write(
        fig1,
        ["chrom", "start", "end", *heights],
        [
            ["chrI", "0", "150", "1", "2", "9", "1", "1", "1"],
            ["chrI", "500", "650", "8", "2", "1", "1", "1", "1"],
            ["chrI", "900", "1050", "1", "1", "1", "1", "1", "7"],
        ],
    )
    # No promoter_gene_id_* or associated_* column: the reader must not need them.
    _write(
        fig2,
        ["chrom_ce10", "start_ce10", "end_ce10", "chrom_ce11", "start_ce11", "end_ce11", "annot"],
        [
            ["chrI", "0", "150", "chrI", "10", "161", "putative_enhancer"],
            ["chrI", "500", "650", "chrI", "510", "661", "coding_promoter"],
            ["chrI", "900", "1050", "chrI", "910", "1061", "putative_enhancer"],
        ],
    )
    _write(
        fig4,
        ["chrom_ce10", "start_ce10", "end_ce10", "devel_is_dynamic"],
        [["chrI", "0", "150", "True"], ["chrI", "500", "650", "True"], ["chrI", "900", "1050", "False"]],
    )
    els, funnel = sd.read_elements(fig1, fig2, fig4)
    assert [(e.start, e.kind, e.tissue) for e in els] == [
        (10, "putative_enhancer", "l2"),
        (510, "coding_promoter", "emb"),
    ]
    assert funnel["rows"] == 3 and funnel["dynamic_by_class"] == {
        "putative_enhancer": 1,
        "coding_promoter": 1,
    }
    assert not set(sd.NEVER_READ) & {"chrom_ce10", "annot", "devel_is_dynamic", *heights}


def test_overlaps_any() -> None:
    rows = [{"chrom": "chrI", "start": s, "end": s + 151} for s in (0, 1_000, 2_000)]
    assert sd.overlaps_any(rows, [("chrI", 150, 300), ("chrI", 1_151, 1_200), ("chrII", 0, 5_000)]) == [
        True,
        False,
        False,
    ]


# --------------------------------------------------------------------------------------------- RNA


def test_read_gene_exons_merges_transcripts_of_protein_coding_genes(tmp_path: Path) -> None:
    import gzip

    lines = [
        "I\tWB\tgene\t100\t900\t.\t+\t.\tID=gene:WBGene1;Name=a-1;biotype=protein_coding;gene_id=WBGene1",
        "I\tWB\tmRNA\t100\t900\t.\t+\t.\tID=transcript:T1.1;Parent=gene:WBGene1;biotype=protein_coding",
        "I\tWB\texon\t100\t200\t.\t+\t.\tParent=transcript:T1.1;Name=e1",
        "I\tWB\texon\t500\t900\t.\t+\t.\tParent=transcript:T1.1;Name=e2",
        "I\tWB\tmRNA\t150\t600\t.\t+\t.\tID=transcript:T1.2;Parent=gene:WBGene1;biotype=protein_coding",
        "I\tWB\texon\t150\t300\t.\t+\t.\tParent=transcript:T1.2;Name=e3",
        "I\tWB\tgene\t2000\t2100\t.\t-\t.\tID=gene:WBGene2;Name=r-1;biotype=ncRNA;gene_id=WBGene2",
        "I\tWB\tncRNA\t2000\t2100\t.\t-\t.\tID=transcript:R1;Parent=gene:WBGene2",
        "I\tWB\texon\t2000\t2100\t.\t-\t.\tParent=transcript:R1",
    ]
    p = tmp_path / "g.gff3.gz"
    with gzip.open(p, "wt") as f:
        f.write("##gff-version 3\n" + "\n".join(lines) + "\n")
    ex = sd.read_gene_exons(p)
    assert ex == {"WBGene1": ("chrI", "+", [(99, 300), (499, 900)])}


def _section(cid: int, items: list[tuple[int, int, float]]) -> bytes:
    raw = struct.pack("<IIIIIBxH", cid, items[0][0], items[-1][1], 0, 0, 1, len(items))
    return raw + b"".join(struct.pack("<IIf", s, e, v) for s, e, v in items)


class _Src:
    def __init__(self, blob: bytes):
        self.blob = blob

    def read(self, off: int, size: int) -> bytes:
        return self.blob[off : off + size]

    def close(self) -> None:
        pass


class _Header:
    uncompress_buf_size = 1 << 16


class _FakeBigWig:
    """A bigWig of zlib-compressed bedGraph sections, enough for chrom_values."""

    def __init__(self, sections: dict[str, list[list[tuple[int, int, float]]]], sizes: dict[str, int]):
        blob, self.items = b"", {}
        self.chroms = {c: (i, sizes[c]) for i, c in enumerate(sizes)}
        for c, secs in sections.items():
            cid = self.chroms[c][0]
            for items in secs:
                data = zlib.compress(_section(cid, items))
                self.items.setdefault(c, []).append(
                    LeafItem(cid, items[0][0], cid, items[-1][1], len(blob), len(data))
                )
                blob += data
        self.src = _Src(blob)
        self.header = _Header()

    def leaf_items(self, chrom: str, intervals: list[tuple[int, int]]) -> list[LeafItem]:
        return list(self.items.get(chrom, []))

    def close(self) -> None:
        pass


def test_chrom_values_fills_bedgraph_runs_and_leaves_gaps_zero() -> None:
    bw = _FakeBigWig({"chrI": [[(0, 3, 1.0), (5, 8, 2.5)], [(8, 10, 4.0)]]}, {"chrI": 12, "chrII": 5})
    v = sd.chrom_values(bw, "chrI")
    assert v.tolist() == [1.0, 1.0, 1.0, 0.0, 0.0, 2.5, 2.5, 2.5, 4.0, 4.0, 0.0, 0.0]
    assert sd.chrom_values(bw, "chrII").tolist() == [0.0] * 5
    assert sd.exon_sums(v, [(0, 2), (7, 9)]) == (2.0 + 2.5 + 4.0, 4)


def test_coverage_reads_each_strand_track(monkeypatch: pytest.MonkeyPatch) -> None:
    tracks = {
        ("emb", "fwd"): _FakeBigWig({"chrI": [[(0, 10, 2.0)]]}, {"chrI": 20}),
        ("emb", "rev"): _FakeBigWig({"chrI": [[(10, 20, 3.0)]]}, {"chrI": 20}),
    }
    monkeypatch.setattr(sd, "BigWig", lambda key: tracks[key])
    exons = {"G1": ("chrI", "+", [(0, 5)]), "G2": ("chrI", "-", [(5, 15)])}
    cov = sd.coverage({k: k for k in tracks}, exons)  # the fake is keyed by its path
    assert cov[("emb", "fwd")] == {"G1": 2.0, "G2": 1.0}
    assert cov[("emb", "rev")] == {"G1": 0.0, "G2": 1.5}


def _cov(sense: float, anti: float) -> tuple[dict, dict]:
    exons = {"P": ("chrI", "+", [(0, 1)]), "M": ("chrI", "-", [(5, 6)])}
    cov = {}
    for s in sd.STAGES:
        cov[(s, "fwd")] = {"P": sense, "M": anti}
        cov[(s, "rev")] = {"P": anti, "M": sense}
    return cov, exons


def test_orientation_check_three_verdicts() -> None:
    assert sd.orientation(*_cov(10.0, 1.0))["verdict"] == "as_labelled"
    assert sd.orientation(*_cov(1.0, 10.0))["verdict"] == "swapped"
    cov, exons = _cov(10.0, 1.0)
    cov[("l3", "fwd")] = {"P": 0.0, "M": 10.0}
    cov[("l3", "rev")] = {"P": 10.0, "M": 0.0}
    assert sd.orientation(cov, exons)["verdict"] == "failed"


def test_expression_is_tpm_on_the_own_strand_and_follows_a_swap() -> None:
    cov, exons = _cov(30.0, 10.0)
    cov[("emb", "fwd")]["P"] = 90.0
    e = sd.expression(cov, exons, "as_labelled", {"P", "M", "absent"})
    assert set(e) == {"P", "M"}
    for i in range(len(sd.STAGES)):
        assert abs(e["P"][i] + e["M"][i] - 1e6) < 1e-6
    assert e["P"][0] == pytest.approx(1e6 * 90 / 120)
    s = sd.expression(cov, exons, "swapped", {"P", "M"})
    assert s["P"][1] == pytest.approx(1e6 * 10 / 20)


def test_concordant_is_stage_ones_rule_on_six_stages() -> None:
    assert sd.concordant((1, 2, 30, 4, 5, 6), "l2")
    assert not sd.concordant((1, 2, 4.9, 4, 0, 0), "l2")  # under 5
    assert not sd.concordant((1, 30, 30, 4, 5, 6), "l2")  # a tie is not highest
    assert not sd.concordant(None, "emb")


def test_stage_control_index_matches_a_brute_force_draw() -> None:
    rng = np.random.default_rng(7)
    genes = [
        simpler.Gene(f"G{i}", f"g{i}", "chrI", p, p + 500, "+" if rng.random() < 0.5 else "-")
        for i, p in enumerate(sorted(rng.integers(0, 2_000_000, 400).tolist()))
    ]
    chroms = {"chrI": simpler.Chromosome(genes)}
    expr = {g.gid: tuple(float(x) for x in rng.integers(0, 20, 6)) for g in genes}
    ctl = sd.ControlIndex(chroms, expr)
    for _ in range(25):
        mid = int(rng.integers(0, 2_000_000))
        g = chroms["chrI"].nearest_tss(mid)
        stage = sd.STAGES[int(rng.integers(0, 6))]
        dist = abs(g.tss - mid)
        b = simpler.distance_bin(dist)
        cands = [h for h in genes if h.gid != g.gid and simpler.distance_bin(abs(h.tss - mid)) == b]
        if not cands:
            continue
        want = sum(sd.concordant(expr[h.gid], stage) for h in cands) / len(cands)
        assert ctl.probability("chrI", stage, mid, g.gid, dist) == pytest.approx(want)


# --------------------------------------------------------------------------------------------- the dose


def test_by_edges_is_the_stage_two_rule() -> None:
    rows = [{"measured": 150, "above": a, "x": a} for a in (0, 5, 12, 20, 28, 40, 60)]
    rows.append({"measured": 0, "above": 0, "x": -1})
    got = sd.by_edges(rows, [0.08, 0.190476])
    assert [(r["x"], r["tier"]) for r in got] == [
        (0, "not_constrained"),
        (5, "not_constrained"),
        (40, "constrained"),
        (60, "constrained"),
    ]
    assert got == st.by_fraction_tertile(rows, [0.08, 0.190476])
    mid = sd.middle_vs_bottom(rows, [0.08, 0.190476])
    assert [(r["x"], r["tier"]) for r in mid] == [
        (0, "not_constrained"),
        (5, "not_constrained"),
        (12, "constrained"),
        (20, "constrained"),
        (28, "constrained"),
    ]


def _population(n: int, seed: int, dose: float, promoter_signal: bool = True):
    """Synthetic genes and elements. Each element's gene is concordant with a base rate set by its
    stage; elements in the top third of constrained fraction add `dose`."""
    rng = np.random.default_rng(seed)
    genes, rows, prom, expr = [], [], [], {}
    for i in range(n):
        pos = 10_000 * i
        g = simpler.Gene(f"E{i}", f"e{i}", "chrI", pos + 300, pos + 2_300, "+")
        genes.append(g)
        stage = sd.STAGES[i % 2]  # emb high base rate, l1 low
        above = int(rng.integers(0, 60))
        top = above / 150 >= 0.2
        base = 0.5 if stage == "emb" else 0.2
        hit = rng.random() < base + (dose if top else 0.0)
        vals = [1.0] * 6
        vals[sd.STAGES.index(stage)] = 50.0 if hit else 0.5
        expr[g.gid] = tuple(vals)
        rows.append(
            {
                "chrom": "chrI",
                "start": pos,
                "end": pos + 151,
                "mid": pos + 75,
                "kind": sd.ELEMENT_TYPE,
                "tissue": stage,
                "domain_gene": g.gid,
                "domain_case": "faces_one",
                "domain_dist": 225,
                "nearest_gene": g.gid if i % 5 else f"P{i}",
                "nearest_dist": 225,
                "measured": 150,
                "above": above,
                "tier": st.tier(150, above),
                "length_measured": 150,
            }
        )
    for i in range(n):
        pos = 10_000 * i + 5_000
        g = simpler.Gene(f"P{i}", f"p{i}", "chrI", pos, pos + 1_000, "+")
        genes.append(g)
        stage = sd.STAGES[i % 6]
        vals = [1.0] * 6
        if promoter_signal or rng.random() < 0.05:
            vals[sd.STAGES.index(stage)] = 50.0
        expr[g.gid] = tuple(vals)
        prom.append(
            {
                "chrom": "chrI",
                "mid": pos - 50,
                "kind": sd.PROMOTER_TYPE,
                "tissue": stage,
                "domain_gene": g.gid,
                "domain_case": "faces_one",
                "domain_dist": 50,
                "nearest_gene": g.gid,
                "nearest_dist": 50,
            }
        )
    chroms = {"chrI": simpler.Chromosome(genes)}
    overlap = [i % 3 == 0 for i in range(n)]
    return chroms, rows, prom, expr, overlap


_OK = {"verdict": "as_labelled", "per_stage": {}}


def test_census_counts_arms_without_any_rna() -> None:
    chroms, rows, prom, _, overlap = _population(1_200, 1, 0.0)
    cen = sd.census(rows, prom, overlap)
    assert cen["elements"] == 1_200 and cen["measured"] == 1_200
    assert cen["dose_arms"]["top"] + cen["dose_arms"]["bottom"] < 1_200
    assert cen["top_in_support"] >= sd.MIN_PER_ARM and not cen["stop_rule_met"]
    assert cen["overlap_stage2_population"] == 400
    assert "c_domain" not in rows[0]
    small = sd.census(rows[:300], prom, overlap[:300])
    assert small["stop_rule_met"]


def test_score_finds_a_planted_dose_and_reads_the_registered_wording() -> None:
    chroms, rows, prom, expr, overlap = _population(3_000, 2, 0.25)
    cen = json.loads(json.dumps(sd.census(rows, prom, overlap)))
    res = sd.score(rows, prom, overlap, chroms, expr, cen, _OK)
    assert res["positive_control"]["separates"]
    assert res["outcome"] == "DOSE_ADDS" and res["falsifier"] == "NOT_FALSIFIED"
    assert res["primary"]["delta_ci95_points"][0] > 0
    w = sd.REGISTRATION["outcome_wording"]
    assert w["DOSE_ADDS"] in res["reading"] and w["always"] in res["reading"]
    assert res["reading"].startswith(w["common"].split("{")[0])
    assert res["secondary"]["outcome"] in ("DOMAIN_RULE_AHEAD", "NO_DIFFERENCE_DETECTED", "NEAREST_AHEAD")


def test_score_without_a_dose_detects_none() -> None:
    chroms, rows, prom, expr, overlap = _population(3_000, 3, 0.0)
    cen = json.loads(json.dumps(sd.census(rows, prom, overlap)))
    res = sd.score(rows, prom, overlap, chroms, expr, cen, _OK)
    assert res["outcome"] == "NO_DIFFERENCE_DETECTED"
    assert sd.REGISTRATION["outcome_wording"]["NO_DIFFERENCE_DETECTED"] in res["reading"]


def test_score_withholds_on_an_unresponsive_readout_and_a_failed_strand() -> None:
    chroms, rows, prom, expr, overlap = _population(1_500, 4, 0.25, promoter_signal=False)
    cen = json.loads(json.dumps(sd.census(rows, prom, overlap)))
    res = sd.score(rows, prom, overlap, chroms, expr, cen, _OK)
    assert res["outcome"] == "UNRESPONSIVE" and "primary" not in res
    assert res["reading"].startswith("UNRESPONSIVE")
    bad = sd.score(rows, prom, overlap, chroms, expr, cen, {"verdict": "failed", "per_stage": {}})
    assert bad["outcome"] == "INSTRUMENT_FAILED" and "positive_control" not in bad


def test_score_stops_when_underpowered() -> None:
    chroms, rows, prom, expr, overlap = _population(300, 5, 0.25)
    cen = json.loads(json.dumps(sd.census(rows, prom, overlap)))
    res = sd.score(rows, prom, overlap, chroms, expr, cen, _OK)
    assert res["outcome"] == "UNDERPOWERED" and "primary" not in res


def test_registration_names_every_part_the_brief_asks_for() -> None:
    r = sd.REGISTRATION
    for key in (
        "hypothesis",
        "population",
        "adaptations",
        "constraint",
        "dose",
        "matching",
        "falsifier",
        "decision",
        "positive_control",
        "secondary",
        "stop_rule",
        "confound",
        "outcome_wording",
    ):
        assert r[key], key
    assert sd.STAGE2_EDGES == (0.08, 0.190476) and sd.STAGE2_CUT == 2.08
    for k in (
        "DOSE_ADDS",
        "NO_DIFFERENCE_DETECTED",
        "DOSE_BELOW",
        "UNDERPOWERED",
        "UNRESPONSIVE",
        "INSTRUMENT_FAILED",
    ):
        assert k in r["outcome_wording"]


# --------------------------------------------------------------------------------------------- local data


LOCAL = (ROOT / sd.FIG1, ROOT / sd.FIG2, ROOT / sd.FIG4, ROOT / simpler.GFF3)


@pytest.mark.skipif(not all(p.exists() for p in LOCAL), reason="git-ignored Jänes 2018 copies absent")
def test_local_population_matches_the_source_counts() -> None:
    els, funnel = sd.read_elements(*LOCAL[:3])
    assert funnel["rows"] == 42_245
    assert funnel["dynamic_by_class"] == {"putative_enhancer": 13_852, "coding_promoter": 10_199}
    reg = ROOT / "data/results/simpler_dose_registration.json"
    if reg.exists():
        assert json.loads(reg.read_text())["element_funnel"] == json.loads(json.dumps(funnel))
