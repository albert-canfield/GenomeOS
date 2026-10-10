# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, the fly (genomeos/attribution/simpler_fly.py, docs/SIMPLER-GENOMES.md section 8).

The tissue grouping, the marker call, the element classes, the endpoint and the matched control are
checked on synthetic inputs. The node is the worm's rule, imported: its cases are tested in
tests/test_simpler_worm.py and are not repeated here.
"""

from __future__ import annotations

import gzip
import random
from pathlib import Path

import pytest

from genomeos.attribution import simpler
from genomeos.attribution import simpler_fly as sf
from genomeos.attribution.simpler import Chromosome, Gene

ROOT = Path(__file__).resolve().parents[1]


def g(gid: str, start: int, end: int, strand: str) -> Gene:
    return Gene(gid, gid, "chrT", start, end, strand)


def test_the_node_is_the_worms_rule_imported_not_copied() -> None:
    assert not hasattr(sf, "domain_rule")
    assert sf.DISTANCE_BINS == simpler.DISTANCE_BINS
    assert (sf.SEED, sf.BOOTSTRAP, sf.CONTROL_DRAWS, sf.MIN_DISCORDANT, sf.ALPHA) == (
        simpler.SEED,
        simpler.BOOTSTRAP,
        simpler.CONTROL_DRAWS,
        simpler.MIN_DISCORDANT,
        simpler.ALPHA,
    )


def test_every_tissue_is_named_by_both_modalities() -> None:
    assert set(sf.ATAC_TISSUE.values()) == set(sf.TISSUES)
    assert set(sf.RNA_TISSUE.values()) == set(sf.TISSUES)
    assert not set(sf.ATAC_TISSUE) & sf.ATAC_EXCLUDED
    assert not set(sf.RNA_TISSUE) & sf.RNA_EXCLUDED
    assert len(set(sf.TISSUES)) == len(sf.TISSUES) == 14


def test_marker_call_keeps_only_single_tissue_peaks() -> None:
    labels = {
        ("10-12", 1): "Neural",
        ("12-14", 2): "Neural prim.",
        ("10-12", 3): "Muscle",
        ("10-12", 4): "Unknown",
    }
    markers = [
        ("chrT-100-300", "10-12", 1, 0.5, 0.001),  # neural in two windows, two labels, one tissue
        ("chrT-100-300", "12-14", 2, 0.4, 0.01),
        ("chrT-500-700", "10-12", 1, 0.5, 0.001),  # neural and muscle
        ("chrT-500-700", "10-12", 3, 0.5, 0.001),
        ("chrT-900-1100", "10-12", 3, 0.5, 0.001),  # muscle, plus an unknown cluster
        ("chrT-900-1100", "10-12", 4, 0.5, 0.001),
        ("chrT-1300-1500", "10-12", 3, 0.5, 0.001),  # muscle; the negative call is ignored
        ("chrT-1300-1500", "10-12", 1, -0.5, 0.001),
        ("chrT-1700-1900", "10-12", 3, 0.5, 0.2),  # not significant
    ]
    calls, funnel = sf.tissue_calls(markers, labels)
    assert calls == {"chrT-100-300": "neural", "chrT-1300-1500": "muscle"}
    assert funnel["peaks_called_in_two_or_more_tissues"] == 1
    assert funnel["peaks_with_a_call_in_an_excluded_label"] == 1
    assert funnel["marker_rows_not_positive_or_not_significant"] == 2


def test_an_ungrouped_label_is_refused() -> None:
    with pytest.raises(ValueError, match="neither grouped nor excluded"):
        sf.tissue_calls([], {("10-12", 1): "Gonad"})


def test_rna_sheets_are_never_opened() -> None:
    for name in sf.NEVER_READ:
        with pytest.raises(ValueError, match="never read"):
            next(iter(sf._sheet(Path("absent.xlsx"), name)))


def test_element_classes() -> None:
    chroms = {"chrT": Chromosome([g("A", 10_000, 20_000, "+"), g("B", 40_000, 50_000, "-")])}
    tss = {"chrT": [10_000, 30_000, 50_000]}  # 30,000: a non-coding transcript's 5' end
    assert sf.classify("chrT", 10_400, chroms, tss) == sf.PROMOTER
    assert sf.classify("chrT", 49_600, chroms, tss) == sf.PROMOTER
    assert sf.classify("chrT", 10_800, chroms, tss) is None  # 800 bp from a 5' end: neither
    assert sf.classify("chrT", 30_500, chroms, tss) is None  # near the non-coding 5' end
    assert sf.classify("chrT", 15_000, chroms, tss) == sf.DISTAL  # intronic, 5 kb from any 5' end
    assert sf.classify("chrU", 15_000, chroms, tss) is None


def test_concordance_needs_detection_and_a_strict_maximum() -> None:
    base = [1.0] * len(sf.TISSUES)
    hi = list(base)
    hi[sf.TISSUES.index("muscle")] = 9.0
    assert sf.concordant(tuple(hi), "muscle")
    low = list(base)
    low[sf.TISSUES.index("muscle")] = 4.9
    assert not sf.concordant(tuple(low), "muscle")
    tie = list(hi)
    tie[sf.TISSUES.index("yolk")] = 9.0
    assert not sf.concordant(tuple(tie), "muscle")
    assert not sf.concordant(None, "muscle")


def test_expression_is_cpm_over_the_grouped_labels(tmp_path: Path) -> None:
    p = tmp_path / "rna.tsv.gz"
    labs = ["brain", "glia", "somatic muscle", "maternal"]
    with gzip.open(p, "wt") as f:
        f.write("gene\t" + "\t".join(labs) + "\n")
        f.write("geneA\t30\t10\t0\t500\n")
        f.write("geneB\t60\t0\t100\t500\n")
        f.write("cr1\t0\t0\t100\t0\n")  # not protein-coding: counts in totals, not in the result
    expr, summary = sf.read_expression(p, {"geneA": "FBgnA", "geneB": "FBgnB"})
    n, m = sf.TISSUES.index("neural"), sf.TISSUES.index("muscle")
    assert summary["tissue_umi_totals"]["neural"] == 100 and summary["tissue_umi_totals"]["muscle"] == 200
    assert expr["FBgnA"][n] == pytest.approx(400_000) and expr["FBgnA"][m] == 0
    assert expr["FBgnB"][n] == pytest.approx(600_000) and expr["FBgnB"][m] == pytest.approx(500_000)
    assert set(expr) == {"FBgnA", "FBgnB"}


def test_an_ungrouped_rna_label_is_refused(tmp_path: Path) -> None:
    p = tmp_path / "rna.tsv.gz"
    with gzip.open(p, "wt") as f:
        f.write("gene\tbrain\tgonad\n" + "geneA\t1\t1\n")
    with pytest.raises(ValueError, match="neither grouped nor excluded"):
        sf.read_expression(p, {})


def test_control_probability_matches_brute_force() -> None:
    rng = random.Random(5)
    genes = []
    for i in range(300):
        s = rng.randrange(0, 2_000_000)
        genes.append(Gene(f"G{i}", f"G{i}", "chrT", s, s + rng.randrange(200, 20_000), rng.choice("+-")))
    chroms = {"chrT": Chromosome(genes)}
    expr = {x.gid: tuple(rng.choice((0.0, 3.0, 6.0, 12.0, 30.0)) for _ in sf.TISSUES) for x in genes}
    idx = sf.ControlIndex(chroms, expr)
    for _ in range(200):
        mid = rng.randrange(0, 2_000_000)
        tissue = rng.choice(sf.TISSUES)
        att = chroms["chrT"].nearest_tss(mid) if rng.random() < 0.5 else rng.choice(genes)
        dist = abs(att.tss - mid)
        b = simpler.distance_bin(dist)
        cands = [x for x in genes if x.gid != att.gid and simpler.distance_bin(abs(x.tss - mid)) == b]
        got = idx.probability("chrT", tissue, mid, att.gid, dist)
        if cands:
            want = sum(sf.concordant(expr[x.gid], tissue) for x in cands) / len(cands)
            assert got == pytest.approx(want)
        else:
            assert got in (0.0, 1.0)


def test_reading_fills_the_registered_wording_only() -> None:
    res = {
        "primary": {
            "n": 1000,
            "domain_concordant": 300,
            "nearest_concordant": 260,
            "domain_rate_pct": 30.0,
            "nearest_rate_pct": 26.0,
            "delta_points": 4.0,
            "delta_ci95_points": [1.0, 7.0],
        },
        "control": {"p_one_sided": 0.001, "expected_rate": 0.2},
        "census": {"discordant": 300},
        "strata": {"intergenic": {"delta_ci95_points": [-1.0, 5.0]}},
    }
    text = sf.reading(res)
    w = sf.REGISTRATION["outcome_wording"]
    for key in ("DOMAIN_RULE_AHEAD", "replication_DOMAIN_RULE_AHEAD", "host_caveat", "tail"):
        assert w[key] in text
    assert "1,000 D. melanogaster" in text and "+4.00 points, 95% +1.00 to +7.00" in text
    res["primary"]["delta_ci95_points"] = [-1.0, 7.0]
    text = sf.reading(res)
    assert w["replication_NO_DIFFERENCE_DETECTED"] in text and w["host_caveat"] not in text


# ----------------------------------------------------------------------- the committed registration

REGISTRATION = ROOT / "data/results/simpler_fly_registration.json"
INVENTORY = ROOT / "data/results/simpler_fly_inventory.json"
LOCAL_HOW = "Calderon 2022 tables and the Ensembl 104 GFF3 live in the git-ignored data/knowledge/simpler_fly"


def test_committed_registration_is_the_module_registration() -> None:
    import json

    import tracked_paths as tp

    tp.must_be_present(REGISTRATION, was="registration not yet written")
    reg = json.loads(REGISTRATION.read_text())
    assert reg["written_before_any_rna_value_was_read"] is True
    assert reg["registration"] == json.loads(json.dumps(sf.REGISTRATION))
    cen = reg["census"]
    assert cen["elements"] == sum(cen["elements_by_tissue"].values())
    assert cen["discordant"] == sum(cen["discordant_by_case"].values())
    assert cen["discordant"] >= sf.MIN_DISCORDANT and not cen["stop_rule_met"]
    assert cen["host_elements"] + cen["intergenic_elements"] == cen["elements"]
    funnel = reg["population_funnel"]["marker_calls"]
    assert funnel["peaks_with_a_positive_marker_call"] == (
        funnel["peaks_with_a_call_in_an_excluded_label"]
        + funnel["peaks_called_in_two_or_more_tissues"]
        + funnel["peaks_called_in_exactly_one_tissue"]
    )
    classes = reg["population_funnel"]["classes"]
    assert classes["distal"] == cen["elements"] and classes["promoter"] == cen["positive_control_promoters"]


def test_inventory_counts_its_downloads_within_budget() -> None:
    import json

    import tracked_paths as tp

    tp.must_be_present(INVENTORY, was="inventory not yet written")
    inv = json.loads(INVENTORY.read_text())
    assert inv["fetched_bytes"] == sum(f["bytes"] for f in inv["fetched"])
    assert inv["lane_download_bytes"] == inv["fetched_bytes"] + inv["scouting_bytes"] <= inv["budget_bytes"]
    assert all(f["licence"] for f in inv["fetched"])
    assert [f["name"] for f in inv["fetched"] if not f["kept"]] == ["GSE190147_scirnaseq_gene_matrix.mtx.gz"]


@pytest.mark.needs_local_data(str(sf.STORE), how=LOCAL_HOW)
def test_census_recomputes_from_the_registered_inputs() -> None:
    import json

    reg = json.loads(REGISTRATION.read_text())
    for key in ("supplementary_tables", "peaks", "gff3", "rna_by_annotation"):
        entry = reg["registered_inputs"][key]
        assert simpler.sha256(ROOT / entry["path"]) == entry["sha256"]
    chroms, tss, _ = sf.read_annotation(ROOT / sf.GFF3)
    calls, funnel = sf.tissue_calls(sf.read_markers(ROOT / sf.SUPP), sf.read_cluster_labels(ROOT / sf.SUPP))
    els, counts = sf.elements(calls, sf.read_peaks(ROOT / sf.PEAKS), chroms, tss)
    assert {"marker_calls": funnel, "classes": counts} == reg["population_funnel"]
    assert sf.census(simpler.attribute(els, chroms)) == reg["census"]
