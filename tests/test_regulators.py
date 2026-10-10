# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1b: the regulator assignment for compiled rules and its registered network-level test.

Synthetic tests of the pure functions in genomeos/attribution/regulators.py, then guarded tests that
read the git-ignored peak clusters and the committed census when they are present."""

from __future__ import annotations

import json
import random
from pathlib import Path

import pytest
import tracked_paths as tp

from genomeos.attribution import regulators as rg

CENSUS = Path("data/results/regulators_census.json")


def rule(cell="K562", action="activates", start=100, end=200, element="E1", gene="G1"):
    return rg.CompiledRule("p", element, gene, action, cell, "predicted", "chr1", start, end)


# --- names and cells ---------------------------------------------------------------------------------


def test_context_label_matches_the_compiled_spelling():
    assert rg.context_label("HeLa-S3") == "HeLa_S3"
    assert rg.context_label("CD14-positive monocyte") == "CD14_positive_monocyte"
    assert rg.context_label("H1-hESC") == "H1"
    assert rg.context_label("K562") == "K562"


def test_motif_names_split_dimers_and_drop_variants():
    assert rg.motif_factor_names("MAX(var.2)") == ["MAX"]
    assert rg.motif_factor_names("Gata1::Tal1") == ["GATA1", "TAL1"]


# --- occupancy ---------------------------------------------------------------------------------------

LINES = [
    "chr1\t50\t150\tGATA1\t500\tK562,HepG2\n",
    "chr1\t160\t170\tTAL1\t400\tK562\n",
    "chr1\t200\t300\tMAX\t300\tHepG2\n",  # touches the end only: half-open, no shared base
    "chr1\t0\t5000\tCTCF\t900\tH1-hESC\n",  # long, starts first: found through the running maximum
    "chr2\t100\t200\tGATA1\t500\tK562\n",
]


def test_overlap_is_half_open_and_finds_long_clusters():
    cl = rg.parse_clusters(LINES)
    occ = rg.element_occupancy(cl, "chr1", 100, 200)
    assert set(occ) == {"GATA1", "TAL1", "CTCF"}
    assert occ["GATA1"] == frozenset({"K562", "HepG2"})
    assert occ["CTCF"] == frozenset({"H1"})
    assert rg.element_occupancy(cl, "chr1", 300, 400) == {"CTCF": frozenset({"H1"})}
    assert rg.element_occupancy(cl, "chr3", 0, 10) == {}


def test_factor_peaks_are_midpoints_in_the_named_cell():
    cl = rg.parse_clusters(LINES)
    assert rg.factor_peaks(cl, "GATA1", "K562") == {"chr1": [100], "chr2": [150]}
    assert rg.factor_peaks(cl, "MAX", "K562") == {}
    assert cl.factors_in()["K562"] == {"GATA1", "TAL1"}


# --- direction is never a default -------------------------------------------------------------------


def test_direction_needs_curated_evidence_that_agrees():
    go = rg.parse_go_directions(
        [
            "!gaf-version: 2.2\n",
            "UniProtKB\tP1\tGATA1\tenables\tGO:0001228\tPMID:1\tIDA\n",
            "UniProtKB\tP2\tREST\tenables\tGO:0001227\tPMID:1\tIDA\n",
            "UniProtKB\tP3\tYY1\tenables\tGO:0001228\tPMID:1\tIDA\n",
            "UniProtKB\tP3\tYY1\tenables\tGO:0001227\tPMID:1\tIDA\n",
            "UniProtKB\tP4\tZZZ\tNOT|enables\tGO:0001228\tPMID:1\tIDA\n",
        ]
    )
    assert rg.factor_class(go, "GATA1") == "activator"
    assert rg.factor_class(go, "REST") == "repressor"
    assert rg.factor_class(go, "YY1") == "activator_and_repressor"
    assert rg.factor_class(go, "ZZZ") == "no_curated_direction"
    assert rg.rule_direction("activator", "activates")[0] == "activates"
    assert rg.rule_direction("repressor", "inhibits")[0] == "inhibits"
    for fclass, action in [
        ("activator", "inhibits"),
        ("repressor", "activates"),
        ("activator_and_repressor", "activates"),
        ("no_curated_direction", "activates"),
    ]:
        assert rg.rule_direction(fclass, action)[0] == "unknown"


def test_candidates_and_kinds():
    go = {"GATA1": {"activator": True, "repressor": False, "codes": set()}}
    occ = {"GATA1": frozenset({"K562"}), "TAL1": frozenset({"K562"}), "MAX": frozenset({"HepG2"})}
    cands = rg.candidates(rule(), occ, frozenset({"GATA1"}), go)
    kinds = {c["factor"]: c["kind"] for c in cands}
    assert kinds == {
        "GATA1": "occupancy_same_cell_with_motif",
        "TAL1": "occupancy_same_cell",
        "MAX": "occupancy_other_cell",
    }
    direction = {c["factor"]: c["direction"] for c in cands}
    assert direction == {"GATA1": "activates", "TAL1": "unknown", "MAX": "unknown"}
    assert rg.rule_kind(cands, 3) == "occupancy_same_cell_with_motif"
    assert rg.rule_kind([], 3) == "motif_only"
    assert rg.rule_kind([], 0) == "unknown"
    # a rule with no cell recorded gains no source anywhere
    assert {c["kind"] for c in rg.candidates(rule(cell=""), occ, frozenset(), go)} == {"occupancy_other_cell"}


def test_tally_counts_a_rule_once_and_its_pairs():
    t = rg.empty_tally()
    go: dict = {}
    occ = {"GATA1": frozenset({"K562"}), "TAL1": frozenset({"K562"})}
    rg.tally(t, rule(), rg.candidates(rule(), occ, frozenset(), go), 5)
    rg.tally(t, rule(element="E2"), [], 0)
    s = rg.finish_tally(t)
    assert s["rules"] == 2 and s["rules_gaining_a_source"] == 1
    assert s["candidate_pairs_in_own_cell"] == 2
    assert s["by_kind"]["occupancy_same_cell"] == 1 and s["by_kind"]["unknown"] == 1
    assert s["direction_of_candidate_pairs"] == {"unknown": 2}


def test_biolang_rule_only_for_a_known_direction():
    r = rule()
    assert rg.biolang_rule("GATA1", r, "unknown", "x") is None
    text = rg.biolang_rule("GATA1", r, "activates", "basis")
    assert text.startswith("rule GATA1 activates G1 { through: E1; when: cell_type = K562;")


def test_compile_spec_names_the_two_engine_blocks():
    assert "strength" in rg.COMPILE_SPEC and "grammar token" in rg.COMPILE_SPEC
    assert "through" in rg.COMPILE_SPEC


# --- ENCODE file choice ------------------------------------------------------------------------------


def efile(acc, ann, rep, analyses=(), date="2020", assembly="GRCh38"):
    return {
        "accession": acc,
        "output_type": "gene quantifications",
        "assembly": assembly,
        "genome_annotation": ann,
        "biological_replicates": [rep],
        "md5sum": "0" * 32,
        "file_size": 10,
        "status": "released",
        "file_format": "tsv",
        "href": f"/files/{acc}/@@download/{acc}.tsv",
        "date_created": date,
        "analyses": list(analyses),
    }


def test_choose_files_takes_the_newest_shared_annotation():
    kd = {
        "default_analysis": "/a/new/",
        "files": [
            efile("K1", "V24", 1),
            efile("K2", "V29", 1, ["/a/new/"]),
            efile("K3", "V29", 2),
            efile("K4", "V19", 1, assembly="hg19"),
        ],
    }
    ctrl = {"default_analysis": None, "files": [efile("C1", "V24", 1), efile("C2", "V24", 2)]}
    got = rg.choose_files(kd, [ctrl])
    assert got["annotation"] == "V24"
    assert [f["accession"] for f in got["knockdown"]] == ["K1"]
    ctrl["files"].append(efile("C3", "V29", 1))
    got = rg.choose_files(kd, [ctrl])
    assert got["annotation"] == "V29" and sorted(f["accession"] for f in got["knockdown"]) == ["K2", "K3"]
    assert rg.choose_files(kd, []) is None


def test_parse_quantification_skips_par_y_and_spikeins():
    lines = [
        "gene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tTPM\tFPKM\n",
        "ENSG00000000003.14\tT\t1\t1\t1\t12.5\t1\n",
        "ENSG00000002586.19_PAR_Y\tT\t1\t1\t1\t9\t1\n",
        "gSpikein_ERCC-00002\tT\t1\t1\t1\t3\t1\n",
    ]
    assert rg.parse_quantification(lines) == {"ENSG00000000003": 12.5}


# --- the registered analysis ------------------------------------------------------------------------


def test_auroc_ties_and_strata():
    assert rg.auroc([1, 2, 3, 4], [False, False, True, True]) == 1.0
    assert rg.auroc([1, 1, 1, 1], [False, True, False, True]) == 0.5
    assert rg.auroc([1, 2], [True, True]) is None
    # within-stratum only: a labelled gene beats the unlabelled one of its own stratum each time
    y = [10, 11, 0, 1]
    assert rg.stratified_auroc(y, [False, True, False, True], [1, 1, 0, 0]) == 1.0
    assert rg.auroc(y, [False, True, False, True]) == 0.75  # unstratified, the 1 loses to the 10


def test_deciles_are_rank_based():
    s = rg.deciles(list(range(100)))
    assert s[0] == 0 and s[99] == 9 and len(set(s)) == 10
    assert rg.deciles([5.0] * 20) == [0] * 20


def test_universe_excludes_factor_neighbours_and_silent_genes():
    tss = {
        "F": ("chr1", 1000),
        "A": ("chr1", 5000),
        "B": ("chr1", 50_000),
        "C": ("chr2", 1000),
        "D": ("chr2", 9),
    }
    ctrl = {"F": 10.0, "A": 10.0, "B": 10.0, "C": 0.5, "D": 3.0}
    assert rg.universe_of(["F", "A", "B", "C", "D", "Z"], ctrl, "F", tss) == ["B", "D"]


def test_baseline_ranks_by_peaks_then_distance():
    tss = {"A": ("chr1", 100_000), "B": ("chr1", 300_000), "C": ("chr1", 900_000)}
    peaks = {"chr1": [100_100, 100_200, 290_000]}
    near = rg.near_peaks(["A", "B", "C"], tss, peaks)
    assert near["A"] == (2, 100) and near["B"] == (1, 10_000)
    assert rg.baseline_arm(["A", "B", "C"], 2, near) == {"A", "B"}


def test_knockdown_ratio_and_response():
    assert rg.knockdown_ratio("F", {"F": 3.0}, {"F": 10.0}) == pytest.approx(0.3)
    assert rg.knockdown_ratio("F", {"F": 0.0}, {"F": 0.5}) is None
    r = rg.response({"A": 3.0}, {"A": 1.0}, ["A"])
    assert r["A"] == pytest.approx(1.0)


def test_derangement_has_no_fixed_point():
    rng = random.Random(1)
    for _ in range(50):
        p = rg.derangement(["a", "b", "c"], rng)
        assert all(x != y for x, y in zip(p, ["a", "b", "c"], strict=True))
    with pytest.raises(ValueError):
        rg.derangement(["a"], rng)


def synthetic_units(n_units: int, signal: float, seed: int = 7) -> list[dict]:
    """Units whose own arm responds by `signal` and whose arms are disjoint between factors."""
    rng = random.Random(seed)
    genes = [f"g{i}" for i in range(400)]
    arms = {f"F{u}": set(genes[u * 15 : u * 15 + 15]) for u in range(n_units)}
    units = []
    for u in range(n_units):
        own = arms[f"F{u}"]
        ctrl = {g: 2.0 + (i % 50) for i, g in enumerate(genes)}
        resp = {g: rng.random() + (signal if g in own else 0.0) for g in genes}
        base = set(rng.sample([g for g in genes if g not in own], len(own)))
        units.append(
            {
                "factor": f"F{u}",
                "cell": "K562",
                "experiments": [rg.experiment_view(genes, ctrl, resp, own, base)],
                "arms": arms,
            }
        )
    return units


def test_analysis_reads_a_planted_factor_specific_signal(monkeypatch):
    monkeypatch.setattr(rg, "PERMUTATIONS", 200)
    monkeypatch.setattr(rg, "BOOTSTRAP_B", 1000)
    out = rg.analyse(synthetic_units(22, signal=0.6))
    assert out["outcome"] == "rules_beat_peaks"
    assert out["control_p"] < 0.05 and out["ci95"][0] > 0
    assert out["reading"].startswith("the genes the compiled rules route")


def test_analysis_without_signal_is_not_factor_specific(monkeypatch):
    monkeypatch.setattr(rg, "PERMUTATIONS", 200)
    monkeypatch.setattr(rg, "BOOTSTRAP_B", 1000)
    out = rg.analyse(synthetic_units(22, signal=0.0))
    assert out["outcome"] == "not_factor_specific"
    assert "decides nothing" in out["reading"]


def test_analysis_stops_below_the_floor_and_sets_aside_lone_cells(monkeypatch):
    monkeypatch.setattr(rg, "PERMUTATIONS", 50)
    monkeypatch.setattr(rg, "BOOTSTRAP_B", 200)
    out = rg.analyse(synthetic_units(5, signal=0.6))
    assert out["outcome"] == "stopped_below_floor" and out["units"] == 5
    units = synthetic_units(21, signal=0.6)
    units[0]["cell"] = "HepG2"  # alone in its cell: cannot be deranged, so it is set aside first
    out = rg.analyse(units)
    assert out["alone_in_cell"] == [{"factor": "F0", "cell": "HepG2"}]
    assert out["units"] == 20


def test_decide_order():
    assert rg.decide(10, 0.1, 0.05, 0.2, 0.001) == "stopped_below_floor"
    assert rg.decide(30, 0.1, 0.05, 0.2, 0.2) == "not_factor_specific"
    assert rg.decide(30, 0.1, 0.05, 0.2, 0.01) == "rules_beat_peaks"
    assert rg.decide(30, -0.1, -0.2, -0.05, 0.01) == "peaks_beat_rules"
    assert rg.decide(30, 0.0, -0.1, 0.1, 0.01) == "no_difference_detected"


def test_every_outcome_has_its_wording():
    assert set(rg.OUTCOMES) == {
        "stopped_below_floor",
        "not_factor_specific",
        "rules_beat_peaks",
        "no_difference_detected",
        "peaks_beat_rules",
    }
    numbers = {"n": 30, "d": 0.01, "lo": -0.02, "hi": 0.04, "p": 0.01, "a": 0.55, "c": 0.52}
    for k in rg.OUTCOMES:
        assert "{" not in rg.reading(k, **numbers)


# --- guarded: the real clusters and the committed census -----------------------------------------------


@pytest.mark.needs_local_data(
    "data/knowledge/regulators/encRegTfbsClusteredWithCells.hg38.bed.gz",
    how="curl -o data/knowledge/regulators/encRegTfbsClusteredWithCells.hg38.bed.gz "
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/encRegTfbsClustered/"
    "encRegTfbsClusteredWithCells.hg38.bed.gz (125 MB; sha256 in regulators.CLUSTER_SOURCE)",
)
def test_real_clusters_name_the_test_cells():
    cl = rg.load_clusters(chroms=("chr21",))
    fac = cl.factors_in()
    assert "K562" in fac and "HepG2" in fac
    assert "GATA1" in fac["K562"]


def test_committed_census_matches_m1_denominator():
    # the census is committed (c267869), so its absence is a broken checkout, not "not run here"
    tp.must_be_present(CENSUS, was="data/results/regulators_census.json is not written yet")
    d = json.loads(CENSUS.read_text())
    t = d["totals"]
    assert t["rules"] == 440_589  # M1's denominator, data/results/runrules_census.json
    assert sum(t["by_kind"].values()) == t["rules"]
    own = t["by_kind"]["occupancy_same_cell_with_motif"] + t["by_kind"]["occupancy_same_cell"]
    assert t["rules_gaining_a_source"] == own
    with_chip = sum(v["rules"] for v in d["per_cell_with_chip"].values())
    without = sum(v["rules"] for v in d["per_cell_without_chip"].values())
    assert with_chip + without == t["rules"]
