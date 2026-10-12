# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2b: the compiler's reader layer (compile.gate_program) and the gated run's scoring (readergate).

Every program here is synthetic and written in the generator's own line format, so the emitter is
tested without the 400 MB generated store; the census and the run read that store.
"""

from __future__ import annotations

import math

import pytest

from genomeos.attribution import bridge
from genomeos.attribution import compile as cp
from genomeos.attribution import context_evidence as ce
from genomeos.attribution import readergate as rg
from genomeos.genome.reader import PeakIndex
from genomeos.lang.parser import BioLangError, parse

#: The tests below reach `scipy`, by one of the two doors this project opens it through:
#: `combine2.solve_z` (`scipy.optimize.brentq`, the competition form's root find) and
#: `runrules._spearman` (`scipy.stats.spearmanr`, which every registered score reads). No bare
#: environment holds it. The reasoning for naming `compose` -- the extra CI's `test` job syncs, of the
#: four that reach scipy -- and what naming one of four costs is written out at tests/test_combine2.py's
#: `needs_scipy`. The marker moves WHERE these tests run, never what they claim.
needs_scipy = pytest.mark.requires_extra("compose")

CELLS = [cid for _, cid in cp.reader_cells()]

PROGRAM = """module human.noncoding.chrT

# a header line the gate must leave alone
gene TGT { symbol: TGT; evidence: curated "GENCODE v50" symbol only, a stub }
gene SUP_1 { symbol: SUP-1; evidence: curated "GENCODE v50" symbol only, a stub }
element E1 {
  class: enhancer
  locus: chrT:100-200
  targets: TGT
  evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect -0.5 log2 fold change on deletion, probability unavailable
}
rule E1 activates TGT { strength: 0.5; when: cell_type = placenta; evidence: predicted "AlphaGenome deletion, placenta" effect -0.5 log2 fold change on deletion, probability unavailable }
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
rule E1_measured activates TGT { strength: 0.2; when: cell_type = K562; evidence: experimental "CRISPRi, silenced in K562"; confidence: 0.90 }
"""  # noqa: E501 - the generator's own line format, which runs past the limit


class FakeReaders:
    """The two methods gate_program reads, over hand-made peaks: E1 open in K562 and HepG2 only, E2
    open nowhere, and no peak file at all for testis (so every element there is unknown)."""

    def __init__(self) -> None:
        self.open_e1 = {"K562", "HepG2"}

    def _peaks(self, biosample: str) -> list[tuple[int, int, float]] | None:
        if biosample == "testis":
            return None
        peaks = [(10, 20, 1.0), (9000, 9100, 1.0)]
        if biosample in self.open_e1:
            peaks.append((150, 160, 5.0))
        return peaks

    def span(self, biosample: str, chrom: str):
        p = self._peaks(biosample)
        return None if p is None else (min(s for s, _, _ in p), max(e for _, e, _ in p))

    def index(self, biosample: str, chrom: str):
        p = self._peaks(biosample)
        return None if p is None else PeakIndex(p)


def gated(selection: str = cp.SELECTION_READER_CELLS):
    return cp.gate_program(PROGRAM, "chrT", selection=selection, readers=FakeReaders())


def _rule(prefix: str) -> str:
    return next(ln for ln in PROGRAM.splitlines() if ln.startswith(prefix))


# ---- the emitter ------------------------------------------------------------------------------------


def test_the_reader_vocabulary_maps_not_assessable_to_unknown_never_closed():
    assert cp.READER_STATE[ce.STATE_NOT_ASSESSABLE] == "unknown"
    assert cp.READER_STATE[ce.STATE_OPEN] == "open"
    assert cp.READER_STATE[ce.STATE_NOT_OPEN] == "closed"
    assert set(cp.READER_STATE) == set(ce.STATES)
    assert [b for b, _ in cp.reader_cells()] == list(ce.READER_TERMS)
    assert len(CELLS) == 13 and "IMR_90" in CELLS and "CD14_positive_monocyte" in CELLS


def test_a_predicted_rule_is_offered_in_the_reader_cells_and_acts_through_its_element():
    line = _rule("rule E1 ")
    out = cp.gate_rule_line(line, cp.SELECTION_READER_CELLS, CELLS)
    assert f"when: cell_type = {'|'.join(CELLS)}; through: E1; " in out
    assert 'evidence: predicted "AlphaGenome deletion, placenta"' in out  # the compiled cell is kept
    kept = cp.gate_rule_line(line, cp.SELECTION_COMPILED_CELL, CELLS)
    assert "when: cell_type = placenta; through: E1; " in kept
    with pytest.raises(ValueError, match="twice"):
        cp.gate_rule_line(out, cp.SELECTION_READER_CELLS, CELLS)
    with pytest.raises(ValueError, match="selection"):
        cp.gate_rule_line(line, "everywhere", CELLS)
    assert cp.gate_rule_line("gene X { }", cp.SELECTION_READER_CELLS, CELLS) == "gene X { }"


def test_a_measured_rule_keeps_the_cell_its_screen_silenced_in():
    line = _rule("rule E1_measured ")
    for sel in cp.SELECTIONS:
        out = cp.gate_rule_line(line, sel, CELLS)
        assert "when: cell_type = K562; through: E1_measured; " in out


def test_the_gated_program_parses_and_states_every_element_in_every_cell():
    text, info = gated()
    m = parse(text, "t")
    assert text.startswith("module human.noncoding.chrT.reader\n")
    assert info["rules_gated"] == {"predicted": 2, "experimental": 1}
    assert info["elements"] == 3
    assert {c.id for c in m.cell_types()} == set(CELLS)
    assert {rd.cell_type for rd in m.readers} == set(CELLS)
    assert all(r.through == r.source for r in m.rules)
    st = info["states_by_cell"]
    assert st["K562"] == {"open": 2, "closed": 1}  # E1 and E1_measured share a locus
    assert st["GM12878"] == {"closed": 3}
    assert st["testis"] == {"unknown": 3}  # no peak file: unknown, never closed
    # every non-rule line of the program survives byte for byte
    body = [ln for ln in PROGRAM.splitlines()[1:] if not ln.startswith("rule ")]
    assert all(ln in text.splitlines() for ln in body)


def test_the_engine_runs_a_rule_only_where_the_reader_states_its_element_open():
    m = parse(gated()[0], "t")
    gate = {c: {r.source: s for r, s, _ in m.reader_gate({"cell_type": c})} for c in CELLS}
    assert gate["K562"] == {"E1": "open", "E2": "closed", "E1_measured": "open"}
    assert gate["HepG2"] == {"E1": "open", "E2": "closed"}  # the measured rule names K562 only
    assert gate["testis"] == {"E1": "unknown", "E2": "unknown"}
    held = [r for r in m.active_rules({"cell_type": "GM12878"}) if r.strength == 0.0]
    assert {r.source for r in held} == {"E1", "E2"}
    # in the compiled cell, which no reader describes, nothing is selected
    assert m.active_rules({"cell_type": "placenta"}) == []


def test_the_compiled_cell_selection_can_only_hold_a_rule_in_its_own_cell():
    m = parse(gated(cp.SELECTION_COMPILED_CELL)[0], "t")
    assert {r.source for r in m.active_rules({"cell_type": "HepG2"})} == set()
    k562 = {r.source: r.strength for r in m.active_rules({"cell_type": "K562"})}
    assert k562 == {"E2": 0.0, "E1_measured": 0.2}  # E2 closed in K562: held, not dropped
    placenta = m.active_rules({"cell_type": "placenta"})
    assert [(r.source, r.strength) for r in placenta] == [("E1", 0.0)]  # unknown there: held


def test_the_gate_refuses_what_it_cannot_gate():
    text, _ = gated()
    with pytest.raises(ValueError, match="already carries a reader"):
        cp.gate_program(text, "chrT", readers=FakeReaders())
    stray = PROGRAM.replace("rule E2 inhibits", "rule E9 inhibits")
    with pytest.raises(ValueError, match="does not declare"):
        cp.gate_program(stray, "chrT", readers=FakeReaders())
    odd = PROGRAM + 'rule E1 activates TGT { strength: 0.5; evidence: predicted "x" }\n'
    with pytest.raises(ValueError, match="were gated"):
        cp.gate_program(odd, "chrT", readers=FakeReaders())
    with pytest.raises(ValueError, match="module line"):
        cp.gate_program("# no module\n" + PROGRAM, "chrT", readers=FakeReaders())


def test_the_reader_blocks_carry_their_source_and_claim_no_function():
    text, _ = gated()
    readers = [ln for ln in text.splitlines() if ln.startswith("reader ")]
    assert len(readers) == 13
    assert all('evidence: experimental "' in ln for ln in readers)
    assert "poised" not in "".join(readers)  # no element-level poised call exists
    assert "accessibility, not function" in text
    assert "silencer" not in "".join(readers).lower()


def test_element_reading_is_the_context_evidence_call():
    fr = FakeReaders()
    assert cp.element_reading("K562", "chrT", 100, 200, fr) == ce.STATE_OPEN
    assert cp.element_reading("GM12878", "chrT", 100, 200, fr) == ce.STATE_NOT_OPEN
    assert cp.element_reading("GM12878", "chrT", 20000, 20100, fr) == ce.STATE_NOT_ASSESSABLE  # beyond span
    assert cp.element_reading("testis", "chrT", 100, 200, fr) == ce.STATE_NOT_ASSESSABLE


# ---- the gated run ----------------------------------------------------------------------------------


def _run(gene: str = "TGT"):
    text, _ = gated()
    # the measured rule would give TGT two rules; the single-rule population never holds such a gene
    text = "\n".join(ln for ln in text.splitlines() if not ln.startswith("rule E1_measured "))
    m = parse(text, "t")
    return rg.run_gene(rg.gated_gene_module(m, gene), gene, rates={}, borrowed=10.0)


def test_the_gated_run_gives_two_levels_and_they_are_the_clamped_ones():
    out = _run()
    assert out is not None
    assert out["states"]["K562"] == "open" and out["states"]["GM12878"] == "closed"
    assert out["states"]["testis"] == "unknown"
    assert out["gate_mismatch"] == []
    assert abs(out["m"] - out["closed_form"]) < rg.CLOSED_FORM_TOLERANCE
    assert out["closed_form"] == pytest.approx(0.5)  # an activating rule: higher where open
    assert out["levels"]["K562"] == pytest.approx(out["L_open"])
    assert out["levels"]["GM12878"] == pytest.approx(out["L_held"])
    assert out["levels"]["testis"] == pytest.approx(out["L_held"])
    assert out["action"] == "activates"


def test_an_inhibiting_rule_is_lower_where_its_element_is_open():
    out = _run("SUP_1")
    assert out is not None and out["action"] == "inhibits"
    assert out["m"] == pytest.approx(-0.3)
    assert set(out["states"].values()) == {"closed", "unknown"}
    assert out["gate_mismatch"] == []
    assert out["symbol"] == "SUP-1"


def test_the_bridge_drops_through_so_the_run_puts_it_back():
    m = parse(gated()[0], "t")
    small = rg.gated_gene_module(m, "SUP_1")
    b = bridge.parameterize(small, {"cell_type": "GM12878"}, build=True, rates={}, borrowed_rate=10.0)
    assert [r.through for r in b.module.rules if r.target == "SUP_1"] == [""]


def test_single_rule_genes_counts_rule_lines_over_programs():
    two = (
        "module x\n"
        'rule E7 activates SUP_1 { strength: 0.1; when: cell_type = K562; evidence: predicted "y" }\n'
    )
    singles, why = rg.single_rule_genes({"chrT": PROGRAM})
    assert singles == {"SUP_1": ("chrT", "E2")}  # TGT has a predicted and a measured rule
    assert why["several_rules"] == 1
    singles, why = rg.single_rule_genes({"chrT": PROGRAM, "chrU": two})
    assert singles == {} and why["genes_with_any_rule"] == 2
    assert rg.compiled_cells(PROGRAM)[("E1", "TGT")] == "placenta"


# ---- the score ---------------------------------------------------------------------------------------


def _gene(i: int, m: float, opened: list[bool], cell: str = "placenta") -> dict:
    return {
        "gene": f"G{i}", "symbol": f"S{i}", "element": f"E{i}",
        "action": "activates" if m > 0 else "inhibits", "m": m, "rho": 2 ** (-m),
        "tier": "borrowed_median", "chrom": "chrT", "compiled_cell": cell,
        "states": {c: ("open" if o else "closed") for c, o in zip(CELLS, opened, strict=True)},
    }  # fmt: skip


def _population(n: int, signal: float, seed: int = 1):
    import random

    rnd = random.Random(seed)
    genes, tpm, ids = [], {c: {} for c in CELLS}, {}
    for i in range(n):
        opened = [rnd.random() < 0.4 for _ in CELLS]
        opened[0], opened[1] = True, False
        m = rnd.choice([0.4, -0.4])
        genes.append(_gene(i, m, opened))
        ids[f"S{i}"] = f"ENSG{i}"
        for c, o in zip(CELLS, opened, strict=True):
            y = 3.0 + signal * (m if o else 0.0) + rnd.gauss(0, 0.5)
            tpm[c][f"ENSG{i}"] = 2**y - 1
    return genes, tpm, ids


def test_assemble_centres_within_the_gene_and_applies_the_population_rules():
    genes, tpm, ids = _population(5, 1.0)
    genes.append(_gene(99, 0.4, [True] * 13))  # open everywhere: no cross-cell prediction
    genes.append(_gene(98, 0.4, [True] + [False] * 12))  # no gene_id for S98
    rows, why = rg.assemble(genes, tpm, ids)
    assert len(rows) == 5
    assert why == {"element_open_in_all_or_none_of_the_cells": 1, "symbol_not_on_exactly_one_gene_id": 1}
    for r in rows:
        assert abs(sum(r["D"])) < 1e-9 and abs(sum(r["P_gated"])) < 1e-9
        assert r["P_ungated"] == [0.0] * 13  # compiled cell placenta is not a reader cell


@needs_scipy
def test_a_gated_prediction_that_is_the_truth_beats_and_one_that_is_not_does_not():
    rows, _ = rg.assemble(*_population(300, 1.0))
    s = rg.score(rows)
    assert s["outcome"] == "BEATS_AND_CONTROL", s["reading"]
    assert s["reading"].startswith("beats baseline: the reader-gated run lowered")
    assert s["ungated_as_compiled"]["clause"] == rg.UNGATED_CLAUSES["closer"]
    rows, _ = rg.assemble(*_population(300, 0.0))
    assert rg.score(rows)["outcome"] in ("WORSE", "NO_DIFFERENCE")
    rows, _ = rg.assemble(*_population(300, -1.0))
    assert rg.score(rows)["outcome"] == "WORSE"


def test_a_small_population_is_not_testable():
    rows, _ = rg.assemble(*_population(rg.MIN_GENES - 1, 1.0))
    s = rg.score(rows)
    assert s["outcome"] == "NOT_TESTABLE" and "fewer than the registered 100" in s["reading"]


def test_every_outcome_sentence_fills_from_the_registered_slots_only():
    nums = {
        "d": "+0.1", "lo": "+0.05", "hi": "+0.2", "n": "1,300", "g": "100", "p": "0.001",
        "versus": rg.UNGATED_CLAUSES["same"], "u": "+0.1", "ulo": "-0.1", "uhi": "+0.2",
        "bad": 1, "min": rg.MIN_GENES,
    }  # fmt: skip
    for kind in rg.OUTCOMES:
        text = rg.reading(kind, **nums)
        assert "{" not in text and "}" not in text
    assert rg.decide(0.1, 0.2, 0.5) == "BEATS_NOT_CONTROL"
    assert rg.decide(-0.2, -0.1, 0.0) == "WORSE"
    assert rg.versus(-0.1, 0.1) == rg.UNGATED_CLAUSES["same"]
    assert rg.versus(-0.2, -0.1) == rg.UNGATED_CLAUSES["further"]


def test_the_registered_words_state_the_independence_caveat_and_claim_no_function():
    assert rg.LIMITATIONS[0].startswith("independence:")
    assert "DNase" in rg.LIMITATIONS[0] and "trained on" in rg.LIMITATIONS[0]
    joined = " ".join([*rg.LIMITATIONS, *rg.OUTCOMES.values(), rg.QUESTION]).lower()
    assert "not evidence that it does" in joined
    for word in ("silencer", "validat", "proves", "confirms"):
        assert word not in joined, word
    assert math.isclose(rg.PSEUDOCOUNT, 1.0)
    assert rg.SELECTION == cp.SELECTION_READER_CELLS


def test_a_gated_program_still_cannot_be_located():
    text, _ = gated()
    with pytest.raises(BioLangError):
        parse(text + 'compartment nucleus { evidence: none "x" }\n', "t")
