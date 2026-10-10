# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1c: the measured CRISPRi pairs traced through each compiled stage (genomeos/attribution/signaltrace.py).

Every program here is synthetic, in the compiled generator's own line format, so each stage is checked
on the committed code it calls (compile.gate_program, bridge.parameterize, readergate.run_gene,
regulators.candidates) without the 450 MB generated store; the run reads that store.
"""

from __future__ import annotations

import math
import random
import re

import pytest

from genomeos.attribution import crispri
from genomeos.attribution import readergate as rgate
from genomeos.attribution import regulators as regs
from genomeos.attribution import signaltrace as st
from genomeos.genome.reader import PeakIndex
from genomeos.results import RESULTS_DIR

PROGRAM = """module human.noncoding.chrT

gene TGT { symbol: TGT; evidence: curated "GENCODE v50" symbol only, a stub }
gene KGN { symbol: KGN; evidence: curated "GENCODE v50" symbol only, a stub }
gene TWO { symbol: TWO; evidence: curated "GENCODE v50" symbol only, a stub }
gene SUP_1 { symbol: SUP-1; evidence: curated "GENCODE v50" symbol only, a stub }
element E1 {
  class: enhancer
  locus: chrT:100-200
  targets: TGT
  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect -0.5 log2 fold change on deletion, probability unavailable
}
rule E1 activates TGT { strength: 0.5; when: cell_type = placenta; evidence: predicted "AlphaGenome deletion, placenta" effect -0.5 log2 fold change on deletion, probability unavailable }
element E4 {
  class: enhancer
  locus: chrT:1000-1100
  targets: KGN
  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect -1.0 log2 fold change on deletion, probability unavailable
}
rule E4 activates KGN { strength: 1.0; when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" effect -1.0 log2 fold change on deletion, probability unavailable }
element E5 {
  class: enhancer
  locus: chrT:2000-2100
  targets: TWO
  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect -0.4 log2 fold change on deletion, probability unavailable
}
rule E5 activates TWO { strength: 0.4; when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" effect -0.4 log2 fold change on deletion, probability unavailable }
element E6 {
  class: enhancer
  locus: chrT:2200-2300
  targets: TWO
  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect -0.2 log2 fold change on deletion, probability unavailable
}
rule E6 activates TWO { strength: 0.2; when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" effect -0.2 log2 fold change on deletion, probability unavailable }
element E2 {
  class: unknown
  locus: chrT:5000-5100
  targets: SUP_1
  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect +0.3 log2 fold change on deletion, probability unavailable
}
rule E2 inhibits SUP_1 { strength: 0.3; when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" effect +0.3 log2 fold change on deletion, probability unavailable }
element E1_measured {
  class: enhancer
  locus: chrT:100-200
  targets: TGT
  evidence: experimental "CRISPRi" links TGT in K562 from 1 pairs
  confidence: 0.90
}
rule E1_measured activates TGT { strength: 0.2; when: cell_type = K562; evidence: experimental "CRISPRi, silenced in K562, held-out split, evaluation only, never a feature"; confidence: 0.90 }
"""  # noqa: E501 - the generator's own line format, which runs past the limit


class FakeReaders:
    """E1 and E4 open in K562 (E1 also in HepG2); E2, E5 and E6 open nowhere; no peak file for testis."""

    def _peaks(self, biosample: str) -> list[tuple[int, int, float]] | None:
        if biosample == "testis":
            return None
        peaks = [(10, 20, 1.0), (9000, 9100, 1.0)]
        if biosample in ("K562", "HepG2"):
            peaks.append((150, 160, 5.0))
        if biosample == "K562":
            peaks.append((1040, 1060, 5.0))
        return peaks

    def span(self, biosample: str, chrom: str):
        p = self._peaks(biosample)
        return None if p is None else (min(s for s, _, _ in p), max(e for _, e, _ in p))

    def index(self, biosample: str, chrom: str):
        p = self._peaks(biosample)
        return None if p is None else PeakIndex(p)


CLUSTERS = [
    "chrT\t150\t160\tGATA1\t1000\tK562,HepG2",  # occupies E1 in K562
    "chrT\t1050\t1060\tCTCF\t1000\tHepG2",  # occupies E4, but not in K562
]


def pair(
    start: int, end: int, gene: str, top: float = 0.0, drop: float = 0.0, chrom: str = "chrT"
) -> crispri.Pair:
    p = crispri.Pair(chrom, start, end, gene, "K562", "test", 5000.0, 1.0, 1.0, False)
    p.features = {"top_target": top, "deletion_drop": drop}
    return p


PAIRS = {
    "tgt": pair(150, 170, "TGT", 1.0, 0.4),  # E1: compiled in placenta, open in K562, GATA1 in K562
    "kgn": pair(1050, 1070, "KGN", 1.0, 0.9),  # E4: compiled in K562, open, no K562 factor
    "two": pair(2050, 2060, "TWO", 1.0, 0.3),  # E5 and E6 both reach TWO: not identifiable
    "sup": pair(5050, 5060, "SUP-1", 0.0, 0.0),  # E2: inhibits, closed in K562
    "other": pair(150, 170, "OTHER", 0.0, 0.2),  # E1 names another gene
    "none": pair(7000, 7100, "TGT", 0.0, 0.0),  # no compiled element
}


@pytest.fixture(scope="module")
def look() -> st.ChromLookup:
    return st.chromosome_lookup(
        "chrT", PROGRAM, PAIRS.values(), {}, 10.0, None, None,
        regs.parse_clusters(CLUSTERS, ["chrT"]), FakeReaders(), RESULTS_DIR,
    )  # fmt: skip


@pytest.fixture(scope="module")
def traced(look) -> dict[str, tuple[dict[str, float], dict[str, str]]]:
    return {k: st.trace_pair(p, look) for k, p in PAIRS.items()}


# ---- the compiled programme read without its experimental rules ------------------------------------------


def test_the_experimental_rule_lines_are_removed_and_counted_against_rule_targets():
    text, removed = st.predicted_text(PROGRAM)
    assert removed == 1
    assert "E1_measured activates" not in text and "rule E1 activates TGT" in text
    assert [k for _, _, k in rgate.rule_targets(text)] == ["predicted"] * 5
    # a line the two readers disagree on is refused, and so is a line rule_targets cannot read
    two_kinds = PROGRAM.replace(
        'evidence: predicted "AlphaGenome deletion, K562" effect +0.3',
        'evidence: predicted "AlphaGenome, evidence: experimental "x" effect +0.3',
    )
    with pytest.raises(ValueError, match="rule_targets"):
        st.predicted_text(two_kinds)
    with pytest.raises(ValueError, match="cannot read"):
        st.predicted_text(
            PROGRAM.replace('evidence: experimental "CRISPRi, silenced', 'evidence: experimental  "x')
        )


def test_the_experimental_layer_is_recorded_for_the_descriptive_and_read_by_no_stage(look, traced):
    assert look.measured == [(100, 200, "TGT", "K562")]
    assert look.measured_names(150, 170, "TGT") and not look.measured_names(150, 170, "KGN")
    # the bridge never saw the measured rule: TGT's one mechanism is the predicted E1, read from its note
    assert look.bridge_reader[("E1", "TGT")] == pytest.approx(2**-0.5)
    assert look.counts["experimental_rule_lines_removed"] == 1


def test_the_overlap_is_half_open_as_the_deletion_table_searches(look):
    assert look.overlapping(150, 170) == ["E1"]
    assert look.overlapping(200, 300) == []  # end is exclusive
    assert look.overlapping(99, 100) == []
    assert look.overlapping(99, 101) == ["E1"]
    assert look.overlapping(2050, 2250) == ["E5", "E6"]


# ---- each stage on the committed code --------------------------------------------------------------------


def test_target_choice_keeps_a_pair_only_when_an_overlapping_rule_names_the_measured_gene(traced):
    assert traced["tgt"][0]["s2_target_kept"] == 1.0 and traced["tgt"][0]["s2_target_drop"] == 0.4
    assert traced["sup"][0]["s2_target_kept"] == 1.0  # symbol SUP-1, gene id SUP_1
    assert traced["other"][0]["s2_target_kept"] == 0.0 and traced["other"][0]["s2_target_drop"] == 0.0
    assert traced["other"][1]["s2_target"] == "the overlapping elements name another gene"
    assert traced["none"][1]["s2_target"] == "no compiled element overlaps"


def test_the_compiled_value_is_the_rules_strength_and_zero_for_an_inhibiting_rule(traced):
    assert traced["tgt"][0]["s3_compiled_value_drop"] == 0.5
    assert traced["kgn"][0]["s3_compiled_value_drop"] == 1.0
    assert (
        traced["sup"][0]["s3_compiled_value_kept"] == 1.0
        and traced["sup"][0]["s3_compiled_value_drop"] == 0.0
    )


def test_the_bridge_under_reader_cells_offers_every_rule_in_k562_and_refuses_two_mechanisms(look, traced):
    assert set(look.bridge_reader) == {("E1", "TGT"), ("E4", "KGN"), ("E2", "SUP_1")}
    assert look.why_reader[("E5", "TWO")] == "not_identifiable"
    assert traced["two"][0]["s4_bridge_kept"] == 0.0 and traced["two"][1]["s4_bridge"] == "not_identifiable"
    assert traced["kgn"][0]["s4_bridge_drop"] == pytest.approx(1.0)
    assert traced["sup"][0]["s4_bridge_kept"] == 1.0 and traced["sup"][0]["s4_bridge_drop"] == 0.0


def test_the_bridge_as_compiled_runs_a_rule_only_in_its_compiled_cell(look, traced):
    assert set(look.bridge_compiled) == {("E4", "KGN"), ("E2", "SUP_1")}
    assert look.why_compiled[("E1", "TGT")] == "rule_not_active_in_K562"
    assert look.why_compiled[("E5", "TWO")] == "not_identifiable"
    assert traced["tgt"][0]["s4a_bridge_as_compiled_kept"] == 0.0
    assert traced["kgn"][0]["s4a_bridge_as_compiled_drop"] == pytest.approx(1.0)


def test_the_reader_gate_keeps_only_an_element_open_in_k562(look, traced):
    # read for the elements a traced pair reaches; E6 reaches no pair here
    assert look.reader == {"E1": "open", "E4": "open", "E5": "closed", "E2": "closed"}
    assert traced["tgt"][0]["s5_reader_kept"] == 1.0
    assert traced["sup"][0]["s5_reader_kept"] == 0.0 and traced["sup"][1]["s5_reader"] == "reader_closed"


def test_the_regulator_stage_needs_a_factor_cluster_naming_k562(look, traced):
    assert look.factors["E1"] == ["GATA1"] and look.factors["E4"] == []
    assert traced["tgt"][0]["s6_regulator_kept"] == 1.0
    assert traced["kgn"][0]["s6_regulator_kept"] == 0.0
    assert traced["kgn"][1]["s6_regulator"] == "no factor cluster names K562 at the element"


def test_the_run_gives_the_closed_form_where_open_and_nothing_where_held(look, traced):
    assert look.runs[("E1", "TGT")]["change"] == pytest.approx(0.5, abs=st.CLOSED_FORM_TOLERANCE)
    assert look.runs[("E4", "KGN")]["change"] == pytest.approx(1.0, abs=st.CLOSED_FORM_TOLERANCE)
    assert look.runs[("E2", "SUP_1")]["state"] == "closed"
    assert look.runs[("E2", "SUP_1")]["change"] == pytest.approx(0.0, abs=st.CLOSED_FORM_TOLERANCE)
    for k in PAIRS:
        f = traced[k][0]
        assert f["s7_run_kept"] == f["s5_reader_kept"]
        assert f["s7_run_drop"] == pytest.approx(f["s5_reader_drop"], abs=st.CLOSED_FORM_TOLERANCE)


def test_the_harness_passes_a_sound_trace_and_names_each_defect(look, traced):
    rows = []
    for k, p in PAIRS.items():
        p.features.update(traced[k][0])
        rows.append(p)
    assert st.harness({"chrT": look}, rows) == []
    bad = st.ChromLookup("chrU", reader={"E9": "open"})
    bad.runs[("E9", "G")] = {"rho": 0.5, "state": "open", "change": 0.9, "gate_mismatch": []}
    assert any("closed form" in x for x in st.harness({"chrU": bad}, []))
    stray = pair(1, 2, "X")
    stray.features.update({f"{s}_{c}": 0.0 for s in st.STAGES for c in ("kept", "drop")})
    stray.features["s5_reader_kept"] = 1.0
    assert st.harness({}, [stray]) == ["s5_reader keeps 1 pairs that s4_bridge dropped"]


# ---- the readings ----------------------------------------------------------------------------------------


def test_decide_follows_the_registered_rule():
    assert st.decide(0, [0.1, 0.2], [0.0, 0.1]) == "NONE_KEPT"
    assert st.decide(5, [-0.2, -0.01], [-0.3, -0.1]) == "WORSE"
    assert st.decide(5, [0.02, 0.2], [-0.1, -0.01]) == "WEAKENED"
    assert st.decide(5, [0.02, 0.2], [-0.1, 0.01]) == "CARRIES"
    assert st.decide(5, [-0.01, 0.2], [-0.2, -0.01]) == "LOST"
    assert st.decide(5, [-0.01, 0.2], [-0.2, 0.01]) == "UNDECIDED"
    assert st.decide(5, [0.0, 0.2], [-0.2, 0.01]) == "UNDECIDED"  # zero is not above zero


def _row(outcome: str) -> dict:
    return {"outcome": outcome, "reading": f"<{outcome}>"}


def test_the_chain_stops_at_the_first_stage_that_keeps_no_gain_and_names_an_earlier_weakening():
    rows = {s: _row("CARRIES") for s in st.STAGES}
    assert st.chain_verdict(rows)[0] == "CARRIED_TO_THE_RUN"
    rows["s3_compiled_value"] = _row("WEAKENED")
    kind, text, at = st.chain_verdict(rows)
    assert (kind, at) == ("WEAKENED_NOT_LOST", "s3_compiled_value")
    rows["s4_bridge"] = _row("NONE_KEPT")
    rows["s5_reader"] = _row("UNDECIDED")
    kind, text, at = st.chain_verdict(rows)
    assert (kind, at) == ("LOST_AT", "s4_bridge")
    assert text.endswith("; it was first weakened at s3_compiled_value (the value the compiled rule states)")
    rows["s4_bridge"] = _row("UNDECIDED")
    assert st.chain_verdict(rows)[:3:2] == ("UNDECIDED_AT", "s4_bridge")
    rows = {s: _row("CARRIES") for s in st.STAGES}
    rows["s6_regulator"] = _row("LOST")  # a branch never stops the chain
    rows["s4a_bridge_as_compiled"] = _row("NONE_KEPT")
    assert st.chain_verdict(rows)[0] == "CARRIED_TO_THE_RUN"


def test_every_outcome_sentence_fills_from_the_registered_slots_only():
    slots = set(st.stage_numbers({
        "gain_over_activity_distance": {"gain": 0.1, "ci95": [0.0, 0.2]},
        "against_anchor": {"gain": -0.1, "ci95": [-0.2, 0.0]},
        "population": {"pairs": 1918, "kept": 3, "kept_regulated": 1},
        "fit_rows": {"pairs": 9000, "kept": 7},
    }))  # fmt: skip
    for kind, text in st.OUTCOMES.items():
        fields = set(re.findall(r"\{(\w+)\}", text))
        assert fields <= (slots | {"what"}), kind
    assert set(st.OUTCOMES) == {"HARNESS_DEFECT", "ANCHOR", *st.LOSS, "CARRIES", "WEAKENED", "UNDECIDED"}


def test_the_registered_words_keep_their_caveats():
    lim = " ".join(st.LIMITATIONS)
    assert "accessibility and occupancy, not function" in lim
    assert "no stage is a test of a mechanism" in lim
    assert "held-out split, evaluation only, never a feature" in st.PREDICTED_ONLY
    assert st.CHAIN[0] == "s1_anchor" and st.CHAIN[-1] == "s7_run"
    assert all(st.STAGES[s]["follows"] in st.STAGES for s in st.BRANCHES)
    assert st.ANCHOR["gain"] == 0.1361 and st.ANCHOR["ci95"] == [0.0749, 0.2258]


def _population(n_chroms: int = 12, per: int = 12, seed: int = 3) -> list[crispri.Pair]:
    rng = random.Random(seed)
    out = []
    for c in range(n_chroms):
        for i in range(per):
            y = rng.random() < 0.25
            p = crispri.Pair(f"chr{c + 1}", i * 10, i * 10 + 5, f"G{c}_{i}", "K562", "t", 1e4, 1.0, 1.0, y)
            signal = (1.0 if y else 0.0) + rng.gauss(0, 0.6)
            p.features = {
                "log_distance": math.log(1000 + rng.random() * 1e5),
                "log_activity": rng.random(),
                "activity_over_distance": rng.random() - 5,
                "top_target": 1.0,
                "deletion_drop": max(0.0, signal),
            }
            p.covered = True
            out.append(p)
    return out


def test_score_reads_a_stage_that_is_the_anchor_as_the_anchor_and_an_empty_stage_as_none_kept():
    fit, pop = _population(seed=1), _population(seed=2)
    for p in fit + pop:
        f = p.features
        for s in st.STAGES:
            f[f"{s}_kept"], f[f"{s}_drop"] = 0.0, 0.0
        for s in ("s1_anchor", "s2_target", "s3_compiled_value", "s4a_bridge_as_compiled"):
            f[f"{s}_kept"], f[f"{s}_drop"] = f["top_target"], f["deletion_drop"]
    anchor = st.anchor_reading(fit, pop)
    assert anchor["recomputed"]["pairs"] == len(pop) and not anchor["reproduced"]
    out = st.score(fit, pop, anchor)
    s = out["stages"]
    assert s["s2_target"]["against_anchor"]["gain"] == 0.0
    assert s["s2_target"]["gain_over_activity_distance"] == s["s1_anchor"]["gain_over_activity_distance"]
    assert s["s4_bridge"]["outcome"] == "NONE_KEPT" and s["s4_bridge"]["population"]["kept"] == 0
    assert (out["outcome"], out["stage_where_it_stops"]) == ("LOST_AT", "s4_bridge")
    assert out["reading"].startswith("the measured signal is lost at s4_bridge")
    assert s["s1_anchor"]["reading"].startswith("anchor reproduced")
