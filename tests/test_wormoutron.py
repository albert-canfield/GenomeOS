# SPDX-License-Identifier: AGPL-3.0-or-later
"""Area I, the worm's outron caveat: the map parses, the masks mean what they say, the labels hold.

The point of these tests is the exclusion. A mask that silently matched nothing would make the lane
report "the lead survives" while having tested nothing, which is the failure mode the verification
practices of 2026-09-17 were written for: a check that preserves the quantity it watches. So the
tests build small worlds where the right answer is known by hand, and they assert what is excluded
AND what is not.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import simpler
from genomeos.attribution import wormoutron as wo
from genomeos.attribution.simpler import Chromosome, Gene

ROOT = Path(__file__).resolve().parents[1]


def _chrom(*genes: Gene) -> dict[str, Chromosome]:
    return {"chrI": Chromosome(list(genes))}


def _tss(start: int, end: int, strand: str, chrom: str = "chrI") -> wo.Interval:
    return wo.Interval(chrom, start, end, strand, "t")


# ------------------------------------------------------------------------------------- the map file


def test_the_map_parses_into_the_four_feature_classes(tmp_path):
    p = tmp_path / "m.gff3"
    p.write_text(
        "##gff-version 3\n"
        "I\toperon\toperon\t100\t900\t.\t+\t.\tName=CEOP1;genes=WBG1,WBG2\n"
        "I\tdicistronic_mRNA\toperon\t2000\t2500\t.\t-\t.\tName=CEOP2;genes=WBG3,WBG4\n"
        "I\tdeprecated_operon\toperon\t3000\t3500\t.\t+\t.\tName=CEOP3;genes=WBG5\n"
        "I\tSL1\tSL1_acceptor_site\t150\t151\t.\t+\t.\tName=s1\n"
        "I\tSL2\tSL2_acceptor_site\t2100\t2101\t.\t-\t.\tName=s2\n"
        "I\tTSS_region\tTSS_region\t50\t70\t.\t+\t.\tName=r1\n"
    )
    m = wo.read_map(p)
    assert [o.name for o in m["operons"]] == ["CEOP1", "CEOP2"]
    assert m["deprecated_operons_skipped"] == 1
    assert m["operons"][0].members == ("WBG1", "WBG2")
    assert len(m["sl_sites"]) == 2
    assert len(m["tss_regions"]) == 1
    # The WS298 sequences are I..X; everything downstream of here speaks chrI..chrX.
    assert {x.chrom for x in m["tss_regions"] + m["sl_sites"]} == {"chrI"}
    assert m["operons"][0].chrom == "chrI"


def test_the_committed_extract_is_the_six_classes_and_nothing_else():
    """The real map, if the store holds it: the filter kept exactly what the registration names."""
    if not wo.MAP_FILE.exists():
        pytest.skip("the git-ignored WS298 extract is not in this checkout's store")
    m = wo.read_map()
    assert len(m["operons"]) == 1385 + 46
    assert m["deprecated_operons_skipped"] == 111
    assert len(m["sl_sites"]) == 43207 + 13444
    assert len(m["tss_regions"]) == 73499


# --------------------------------------------------------------------------------------- the outron


def test_the_outron_runs_from_the_upstream_tss_region_to_the_annotated_five_prime_end():
    g = Gene("A", "a", "chrI", 1000, 2000, "+")
    out = wo.outrons(_chrom(g), [_tss(700, 720, "+")])
    assert out["A"] == (700, 1000)


def test_the_outron_of_a_minus_strand_gene_lies_above_its_annotated_five_prime_end():
    g = Gene("A", "a", "chrI", 1000, 2000, "-")  # 5' end is 2000
    out = wo.outrons(_chrom(g), [_tss(2200, 2300, "-")])
    assert out["A"] == (2000, 2300)


def test_a_tss_region_on_the_other_strand_gives_no_outron():
    g = Gene("A", "a", "chrI", 1000, 2000, "+")
    assert wo.outrons(_chrom(g), [_tss(700, 720, "-")]) == {}


def test_a_tss_region_beyond_the_cap_gives_no_outron_rather_than_an_outron_of_the_cap():
    g = Gene("A", "a", "chrI", 100_000, 101_000, "+")
    far = _tss(100_000 - wo.OUTRON_MAX_BP - 1, 100_000 - wo.OUTRON_MAX_BP, "+")
    assert wo.outrons(_chrom(g), [far]) == {}
    near = _tss(100_000 - wo.OUTRON_MAX_BP, 100_000 - wo.OUTRON_MAX_BP + 5, "+")
    assert wo.outrons(_chrom(g), [near])["A"] == (100_000 - wo.OUTRON_MAX_BP, 100_000)


def test_the_outron_is_truncated_at_an_intervening_same_strand_five_prime_end():
    upstream = Gene("U", "u", "chrI", 800, 900, "+")  # its 5' end, 800, sits inside the raw span
    g = Gene("A", "a", "chrI", 1000, 2000, "+")
    out = wo.outrons(_chrom(upstream, g), [_tss(700, 720, "+")])
    assert out["A"] == (801, 1000)


def test_the_outron_ignores_a_tss_region_downstream_of_the_five_prime_end():
    g = Gene("A", "a", "chrI", 1000, 2000, "+")
    assert wo.outrons(_chrom(g), [_tss(1500, 1520, "+")]) == {}


# ------------------------------------------------------------------- trans-splice acceptor evidence


def test_sl_evidence_is_the_acceptor_site_at_the_annotated_five_prime_end():
    g = Gene("A", "a", "chrI", 1000, 2000, "+")
    h = Gene("B", "b", "chrI", 5000, 6000, "+")
    sl = [wo.Interval("chrI", 1000 + wo.SL_WINDOW_BP, 1000 + wo.SL_WINDOW_BP + 1, "+", "s")]
    assert wo.sl_spliced_genes(_chrom(g, h), sl) == {"A"}
    far = [wo.Interval("chrI", 1000 + wo.SL_WINDOW_BP + 1, 1000 + wo.SL_WINDOW_BP + 2, "+", "s")]
    assert wo.sl_spliced_genes(_chrom(g, h), far) == set()


# ---------------------------------------------------------------------------------------- the operon


def test_the_first_operon_member_is_the_five_prime_most_one_on_the_operon_strand():
    a = Gene("A", "a", "chrI", 100, 200, "+")
    b = Gene("B", "b", "chrI", 400, 500, "+")
    chroms = _chrom(a, b)
    plus = wo.Operon("chrI", 100, 500, "+", "CEOP1", ("A", "B"), "operon")
    m = wo.operon_members([plus], chroms)
    assert m["A"]["is_first"] and not m["B"]["is_first"]
    minus = wo.Operon("chrI", 100, 500, "-", "CEOP2", ("A", "B"), "operon")
    m = wo.operon_members([minus], chroms)
    assert m["B"]["is_first"] and not m["A"]["is_first"]


def test_a_member_whose_order_cannot_be_known_counts_as_downstream():
    a = Gene("A", "a", "chrI", 100, 200, "+")
    op = wo.Operon("chrI", 100, 500, "+", "CEOP1", ("A", "NOT_IN_MODELS"), "operon")
    m = wo.operon_members([op], _chrom(a))
    assert m["A"]["order_known"] is False
    assert m["A"]["is_first"] is False  # conservative: it is excluded


def test_the_span_index_answers_containment_at_the_edges():
    idx = wo.SpanIndex([(10, 20), (100, 200), (150, 160)])
    assert all(idx.covers(p) for p in (10, 20, 100, 155, 200))
    assert not any(idx.covers(p) for p in (9, 21, 99, 201))


# ----------------------------------------------------------------------------------------- the masks


def _row(mid: int, dom: str, near: str, case: str = "faces_one") -> dict:
    return {
        "chrom": "chrI",
        "mid": mid,
        "kind": simpler.ELEMENT_TYPE,
        "tissue": "Neurons",
        "domain_gene": dom,
        "domain_case": case,
        "domain_dist": 100,
        "nearest_gene": near,
        "nearest_dist": 100,
    }


def test_the_primary_mask_fires_on_an_outron_and_spares_an_element_outside_every_mask():
    g = Gene("A", "a", "chrI", 1000, 2000, "+")
    far = Gene("B", "b", "chrI", 50_000, 51_000, "+")
    chroms = _chrom(g, far)
    mapdata = {"operons": [], "sl_sites": [], "tss_regions": [_tss(700, 720, "+")]}
    inside = _row(800, "A", "A")
    outside = _row(40_000, "B", "B")
    wo.flag([inside, outside], chroms, mapdata)
    assert inside["outron_domain"] == 1 and inside["excluded"] == 1
    assert outside["excluded"] == 0, "a clean element must survive, or the mask tests nothing"


def test_the_operon_mask_fires_on_the_span_and_on_a_downstream_member():
    a = Gene("A", "a", "chrI", 1000, 1100, "+")
    b = Gene("B", "b", "chrI", 1400, 1500, "+")
    away = Gene("C", "c", "chrI", 90_000, 91_000, "+")
    chroms = _chrom(a, b, away)
    op = wo.Operon("chrI", 1000, 1500, "+", "CEOP1", ("A", "B"), "operon")
    mapdata = {"operons": [op], "sl_sites": [], "tss_regions": []}
    in_span = _row(1200, "B", "B")
    on_member = _row(80_000, "B", "C")  # outside the span, but attributed to a downstream member
    clean = _row(80_000, "C", "C")
    wo.flag([in_span, on_member, clean], chroms, mapdata)
    assert in_span["operon_span"] == 1 and in_span["excluded"] == 1
    assert on_member["operon_span"] == 0 and on_member["operon_downstream_domain"] == 1
    assert on_member["excluded"] == 1
    assert clean["excluded"] == 0


def test_the_strict_mask_is_at_least_as_wide_as_the_primary_one():
    g = Gene("A", "a", "chrI", 1000, 2000, "+")
    other = Gene("B", "b", "chrI", 4000, 5000, "+")
    chroms = _chrom(g, other)
    mapdata = {
        "operons": [],
        "sl_sites": [],
        "tss_regions": [_tss(700, 720, "+"), _tss(3500, 3520, "+")],
    }
    r = _row(3800, "A", "A")  # in B's outron, in nobody's attributed outron
    wo.flag([r], chroms, mapdata)
    assert r["excluded"] == 0 and r["excluded_strict"] == 1


# ------------------------------------------------------------------------------- census and labels


def test_the_census_counts_only_between_gene_elements_and_states_its_floors():
    rows = [
        {**_row(10, "A", "A"), "excluded": 0, "excluded_strict": 0, "excluded_sl": 0},
        {**_row(20, "A", "B"), "excluded": 1, "excluded_strict": 1, "excluded_sl": 1},
        {**_row(30, "A", "B", case="host"), "excluded": 0, "excluded_strict": 0, "excluded_sl": 0},
    ]
    for r in rows:
        for k in (
            "outron_either",
            "outron_domain",
            "outron_nearest",
            "outron_either_sl",
            "outron_any_gene",
            "operon_span",
            "operon_downstream_domain",
            "operon_downstream_nearest",
            "operon_mask",
        ):
            r.setdefault(k, 0)
    cen = wo.census(rows, {})
    assert cen["between_genes_recomputed"] == 2  # the host element is not this lane's population
    assert cen["matches_stage1"] is False
    assert cen["surviving"] == 1
    assert cen["minimum_surviving"] == 800
    assert cen["minimum_surviving_discordant"] == 100
    assert cen["inconclusive_unscored"] is True


def test_the_survival_labels_are_the_registered_decision():
    lo_stage1 = wo.STAGE1_BETWEEN_GENES["delta_ci95_points"][0]
    assert wo.survival_label(4.0, 1.0, 7.0) == "LEAD_SURVIVES"
    assert wo.survival_label(lo_stage1, 0.5, 6.0) == "LEAD_SURVIVES"
    assert wo.survival_label(lo_stage1 - 0.01, 0.5, 6.0) == "LEAD_ATTENUATED"
    assert wo.survival_label(1.0, -1.0, 3.0) == "LEAD_DOES_NOT_SURVIVE"
    assert wo.survival_label(-2.0, -4.0, -1.0) == "NEAREST_AHEAD_AFTER_EXCLUSION"


def test_an_inconclusive_run_reads_its_floors_and_scores_nothing():
    cen = {"surviving": 10, "surviving_discordant": 2, "inconclusive_unscored": True}
    res = wo.score([], {}, {}, cen)
    assert res["label"] == "INCONCLUSIVE"
    assert res["primary"] is None
    assert "800 and 100" in res["reading"]
    assert "re-analysis" in res["reading"]


# -------------------------------------------------------------------------- the registration itself


def test_the_registration_states_the_prior_exposure_plainly():
    r = wo.REGISTRATION
    assert "Serizay" in r["prior_exposure"]
    assert "ALREADY READ" in r["prior_exposure"]
    assert "NOT A FRESH TEST" in r["prior_exposure"]
    assert "re-analysis" in r["outcome_wording"]["exposure"]


def test_the_registration_carries_stage_one_wording_verbatim():
    for key in ("baseline", "endpoint", "statistic", "matched_control"):
        assert wo.REGISTRATION[key] == simpler.REGISTRATION[key], key
    assert simpler.REGISTRATION["attribution"] in wo.REGISTRATION["attribution"]


def test_the_registration_names_both_halves_of_the_falsifier_and_the_inconclusive_rule():
    f = wo.REGISTRATION["falsifier"]
    assert "crosses 0" in f and "+2.79" in f
    assert "800" in wo.REGISTRATION["inconclusive"]
    assert "100" in wo.REGISTRATION["inconclusive"]


def test_no_reading_template_claims_regulation_or_function():
    assert "not that the element acts on it" in wo.REGISTRATION["outcome_wording"]["tail"]
    assert "Accessibility is not function" in wo.REGISTRATION["what_it_cannot_establish"]
    for text in wo.REGISTRATION["outcome_wording"].values():
        if "regulates the gene" in text:
            assert "never that" in text or "not that" in text


def test_the_lane_spends_nothing_and_asks_no_model():
    text = (ROOT / "genomeos" / "attribution" / "wormoutron.py").read_text()
    for forbidden in ("alphagenome", "requests.post", "ANTHROPIC_API_KEY", "urllib.request"):
        assert forbidden not in text, forbidden


def test_the_committed_registration_matches_the_module_if_it_exists():
    p = ROOT / "data" / "results" / "wormoutron_registration.json"
    if not p.exists():
        pytest.skip("the registration has not been written yet")
    reg = json.loads(p.read_text())
    assert reg["registration"] == json.loads(json.dumps(wo.REGISTRATION))
