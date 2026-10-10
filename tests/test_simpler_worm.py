# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, C. elegans first (genomeos/attribution/simpler.py, docs/SIMPLER-GENOMES.md).

The rule cases are checked on synthetic genes; the committed registration is checked against the
module, so the registered wording cannot drift; the census is recomputed from the local copies when
they are present (they are git-ignored, so only that test skips without them).
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import pytest

from genomeos.attribution import simpler
from genomeos.attribution.simpler import Chromosome, Gene

ROOT = Path(__file__).resolve().parents[1]
REGISTRATION = ROOT / "data/results/simpler_worm_registration.json"
INVENTORY = ROOT / "data/results/simpler_worm_inventory.json"
RESULT = ROOT / "data/results/simpler_worm_result.json"
LOCAL = (ROOT / simpler.TABLE_S2, ROOT / simpler.GFF3)


def g(gid: str, start: int, end: int, strand: str) -> Gene:
    return Gene(gid, gid, "chrT", start, end, strand)


# A + gene, then a divergent pair, then a convergent pair, then a tandem pair.
CHROM = Chromosome(
    [
        g("A", 1_000, 5_000, "+"),
        g("B", 8_000, 9_000, "-"),  # B and C diverge from the gap 9,000-10,000
        g("C", 10_000, 14_000, "+"),
        g("D", 20_000, 21_000, "-"),  # C and D converge on the gap 14,000-20,000
        g("E", 30_000, 31_000, "+"),  # E then F, tandem
        g("F", 36_000, 40_000, "+"),
        g("H", 3_000, 4_000, "-"),  # nested in A
    ]
)


def test_host_wins_even_when_another_tss_is_nearer() -> None:
    gene, case = simpler.domain_rule(CHROM, 4_800)
    assert (gene.gid, case) == ("A", "host")
    assert CHROM.nearest_tss(4_800).gid == "H"  # H's 5' end is at 4,000


def test_overlapping_hosts_take_the_nearer_5prime_end() -> None:
    gene, case = simpler.domain_rule(CHROM, 3_500)
    assert (gene.gid, case) == ("H", "host")


def test_divergent_pair_takes_the_nearer_5prime_end() -> None:
    assert simpler.domain_rule(CHROM, 9_400)[0].gid == "B"
    assert simpler.domain_rule(CHROM, 9_700)[0].gid == "C"
    assert simpler.domain_rule(CHROM, 9_700)[1] == "divergent"


def test_convergent_pair_falls_back_to_the_nearer_5prime_end() -> None:
    gene, case = simpler.domain_rule(CHROM, 15_000)
    assert case == "convergent_fallback"
    assert gene.gid == "C"  # C's 5' end 10,000 is 5 kb away, D's 21,000 is 6 kb


def test_tandem_pair_takes_the_downstream_gene_it_sits_upstream_of() -> None:
    gene, case = simpler.domain_rule(CHROM, 31_500)
    assert (gene.gid, case) == ("F", "faces_one")
    assert CHROM.nearest_tss(31_500).gid == "E"  # E's 5' end at 30,000 is nearer than F's 36,000


def test_one_flank_at_the_chromosome_end() -> None:
    gene, case = simpler.domain_rule(CHROM, 50_000)
    assert (gene.gid, case) == ("F", "one_flank")


def test_concordance_needs_detection_and_a_strict_maximum() -> None:
    assert simpler.concordant((1.0, 9.0, 2.0, 3.0, 4.0), "Neurons")
    assert not simpler.concordant((1.0, 4.9, 2.0, 3.0, 4.0), "Neurons")  # below 5 TPM
    assert not simpler.concordant((9.0, 9.0, 2.0, 3.0, 4.0), "Neurons")  # tied
    assert not simpler.concordant(None, "Neurons")


def test_control_probability_matches_brute_force() -> None:
    rng = random.Random(3)
    genes = []
    for i in range(300):
        s = rng.randrange(0, 2_000_000)
        genes.append(Gene(f"G{i}", f"G{i}", "chrT", s, s + rng.randrange(200, 20_000), rng.choice("+-")))
    chroms = {"chrT": Chromosome(genes)}
    expr = {x.gid: tuple(rng.choice((0.0, 3.0, 6.0, 12.0, 30.0)) for _ in simpler.TISSUES) for x in genes}
    idx = simpler.ControlIndex(chroms, expr)
    for _ in range(200):
        mid = rng.randrange(0, 2_000_000)
        tissue = rng.choice(simpler.TISSUES)
        att = chroms["chrT"].nearest_tss(mid) if rng.random() < 0.5 else rng.choice(genes)
        dist = abs(att.tss - mid)
        b = simpler.distance_bin(dist)
        cands = [x for x in genes if x.gid != att.gid and simpler.distance_bin(abs(x.tss - mid)) == b]
        got = idx.probability("chrT", tissue, mid, att.gid, dist)
        if cands:
            want = sum(simpler.concordant(expr[x.gid], tissue) for x in cands) / len(cands)
            assert got == pytest.approx(want)
        else:
            assert got in (0.0, 1.0)


def test_cluster_bootstrap_is_seeded_and_brackets_the_point() -> None:
    rows = [{"chrom": "chrT", "mid": i * 50_000, "d": (1 if i % 3 == 0 else 0)} for i in range(400)]
    a = simpler.cluster_bootstrap(rows, "d", reps=300)
    b = simpler.cluster_bootstrap(rows, "d", reps=300)
    assert a == b
    assert a["lo"] <= a["point"] <= a["hi"]


def test_comparison_labels() -> None:
    assert simpler.comparison_label(0.1, 2.0) == "DOMAIN_RULE_AHEAD"
    assert simpler.comparison_label(-2.0, -0.1) == "NEAREST_AHEAD"
    assert simpler.comparison_label(-1.0, 1.0) == "NO_DIFFERENCE_DETECTED"


def test_committed_registration_is_the_module_registration() -> None:
    reg = json.loads(REGISTRATION.read_text())
    assert reg["written_before_the_rna_sheet_was_read"] is True
    assert reg["registration"] == json.loads(json.dumps(simpler.REGISTRATION))
    cen = reg["census"]
    assert cen["elements"] == sum(cen["elements_by_tissue"].values())
    assert cen["discordant"] == sum(cen["discordant_by_case"].values())
    assert cen["discordant"] >= simpler.MIN_DISCORDANT
    assert cen["host_elements"] + cen["intergenic_elements"] == cen["elements"]


def test_inventory_counts_its_downloads() -> None:
    inv = json.loads(INVENTORY.read_text())
    fetched = [s for s in inv["sources"] if s["fetched"]]
    assert len(fetched) == 1
    assert inv["lane_download_bytes"] == fetched[0]["bytes"]["bundle_streamed_then_discarded"]
    assert inv["lane_download_bytes"] <= 1_500_000_000


@pytest.mark.skipif(not all(p.exists() for p in LOCAL), reason="git-ignored Table S2 or WBcel235 GFF3 absent")
def test_census_recomputes_from_the_registered_inputs() -> None:
    reg = json.loads(REGISTRATION.read_text())
    for key in ("table_s2", "gff3"):
        assert (
            simpler.sha256(ROOT / reg["registered_inputs"][key]["path"])
            == reg["registered_inputs"][key]["sha256"]
        )
    chroms = simpler.read_genes(ROOT / simpler.GFF3)
    rows = simpler.attribute(simpler.read_elements(ROOT / simpler.TABLE_S2), chroms)
    assert simpler.census(rows) == reg["census"]


def test_result_reading_is_the_registered_wording() -> None:
    res = json.loads(RESULT.read_text())
    reg = json.loads(REGISTRATION.read_text())
    assert res["census"] == reg["census"]
    assert res["registration_commit"]
    assert res["positive_control"]["responsive"] is True
    assert (
        res["outcome"] == res["comparison"] == simpler.comparison_label(*res["primary"]["delta_ci95_points"])
    )
    assert simpler.reading({**res, "strata": res["strata_beside"]}) == res["reading"]
    assert res["alphagenome_requests"] == 0
