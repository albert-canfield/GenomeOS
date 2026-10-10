# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-combine2: the joint fit of multi-element genes under each combination form, the gated run that
must equal it, and the registered scoring (genomeos/attribution/combine2.py).

Every program here is synthetic and written in the generator's own line format, so nothing reads the
400 MB generated store; the census and the run read that store.
"""

from __future__ import annotations

import ast
import json
import math
import random
from pathlib import Path

import pytest
import tracked_paths as tp

from genomeos.attribution import combine2 as c2
from genomeos.attribution import compile as cp
from genomeos.attribution import readergate as rg
from genomeos.genome.reader import PeakIndex
from genomeos.lang.parser import parse
from genomeos.runtime import grn

CELLS = list(c2.CELLS)

_EV = 'evidence: predicted "AlphaGenome deletion, {cell}" effect {x:+.3f} log2 fold change on deletion, probability unavailable'  # noqa: E501
_EL = 'evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" effect {x:+.3f} log2 fold change on deletion, probability unavailable'  # noqa: E501


def _element(eid: str, start: int, gene: str, x: float) -> str:
    return (
        f"element {eid} {{\n  class: enhancer\n  locus: chrT:{start}-{start + 100}\n  targets: {gene}\n"
        f"  {_EL.format(x=x)}\n}}"
    )


def _rule(eid: str, gene: str, x: float, cell: str) -> str:
    act = "activates" if x < 0 else "inhibits"
    ev = _EV.format(x=x, cell=cell)
    return f"rule {eid} {act} {gene} {{ strength: {abs(x):.3f}; when: cell_type = {cell}; {ev} }}"


#: TGT: two activating elements and one inhibiting one, compiled in three different cells; TWO: two
#: activating elements; ONE: a single rule (not a multi-element gene)
ELEMENTS = {
    "E1": (100, "TGT", -0.5, "placenta"),
    "E2": (300, "TGT", -0.2, "K562"),
    "E3": (500, "TGT", +0.1, "liver"),
    "F1": (700, "TWO", -0.3, "K562"),
    "F2": (900, "TWO", -0.1, "HepG2"),
    "G1": (1100, "ONE", -0.4, "K562"),
}
PROGRAM = (
    "\n".join(
        [
            "module human.noncoding.chrT",
            *(
                f'gene {g} {{ symbol: {g}; evidence: curated "GENCODE v50" symbol only, a stub }}'
                for g in ("TGT", "TWO", "ONE")
            ),
            *(_element(e, s, g, x) for e, (s, g, x, _c) in ELEMENTS.items()),
            *(_rule(e, g, x, c) for e, (_s, g, x, c) in ELEMENTS.items()),
        ]
    )
    + "\n"
)
#: biosample -> the elements a narrowPeak overlaps there; testis has no file (every element unknown)
OPEN = {
    "K562": {"E1", "E2", "F1", "G1"},
    "HepG2": {"E1", "E3", "F2"},
    "GM12878": {"E2"},
    "IMR-90": {"E3", "F1", "F2"},
}
LABEL = dict.fromkeys(("source", "species_cell", "tier"), "t")


class FakeReaders:
    def _peaks(self, biosample: str):
        if biosample == "testis":
            return None
        peaks = [(10, 20, 1.0), (9000, 9100, 1.0)]
        for e in OPEN.get(biosample, ()):
            s = ELEMENTS[e][0]
            peaks.append((s + 40, s + 60, 5.0))
        return peaks

    def span(self, biosample: str, chrom: str):
        p = self._peaks(biosample)
        return None if p is None else (min(s for s, _, _ in p), max(e for _, e, _ in p))

    def index(self, biosample: str, chrom: str):
        p = self._peaks(biosample)
        return None if p is None else PeakIndex(sorted(p))


def _module():
    text, _ = cp.gate_program(PROGRAM, "chrT", selection=c2.SELECTION, readers=FakeReaders())
    return parse(text, "t")


def _genes():
    m = _module()
    return m, {gr.gene: gr for gr in c2.gene_rules(m, {"TGT", "TWO"}, "chrT", rg.compiled_cells(PROGRAM))}


def _rho(*xs: float) -> dict[str, float]:
    return {f"e{i}": 2.0**x for i, x in enumerate(xs)}


# ---- the joint fit ---------------------------------------------------------------------------------------


@pytest.mark.parametrize("form", c2.FORMS)
def test_the_fit_reproduces_every_single_removal_and_the_intact_level(form):
    rho = _rho(-0.5, -0.2, -0.05, +0.1, +0.03)
    f = c2.fit(form, rho)
    assert f.ok, f.reason
    everything = dict.fromkeys(rho, True)
    assert c2.level(f, everything) == pytest.approx(1.0, rel=1e-12)
    for el, r in rho.items():
        assert c2.level(f, {**everything, el: False}) == pytest.approx(r, rel=1e-9)
    assert 0 < f.q < 1
    assert all(0 < s <= 1 for s in f.s.values())
    assert max(f.s[el] for el in f.weights) == 1.0  # DECLARED_STRENGTH convention


def test_the_mean_is_linear_in_the_contributions_under_the_hold():
    rho = _rho(-0.5, -0.2, -0.3)
    f = c2.fit("mean", rho)
    d = {el: 1 - r for el, r in rho.items()}
    assert f.q == pytest.approx(sum(d.values()))
    none = dict.fromkeys(rho, False)
    assert c2.level(f, none) == pytest.approx(1 - f.q)
    for el in rho:  # opening one element adds exactly its own contribution, whatever else is closed
        assert c2.level(f, {**none, el: True}) - c2.level(f, none) == pytest.approx(d[el])


def test_competition_occupancy_is_the_declared_theta_and_its_singles_fall_short_of_adding():
    rho = _rho(-0.5, -0.2, -0.3)
    f = c2.fit("competition", rho)
    k = c2.THETA / (1 - c2.THETA)
    assert sum(f.weights.values()) == pytest.approx(k)
    none = dict.fromkeys(rho, False)
    total = 1 - c2.level(f, none)  # the element-driven share
    singles = sum(1 - r for r in rho.values())
    assert singles < total  # sub-additive: removing all lowers more than the singles add up to
    assert f.q > c2.fit("mean", rho).q


def test_competition_tends_to_the_mean_as_engagement_weakens():
    rho = _rho(-0.4, -0.15, +0.05)
    weak = c2.fit("competition", rho, theta=1e-7)
    mean = c2.fit("mean", rho)
    for pattern in ((True, False, False), (False, True, True), (False, False, False), (True, True, False)):
        o = dict(zip(rho, pattern, strict=True))
        assert c2.level(weak, o) == pytest.approx(c2.level(mean, o), rel=1e-5)


def test_what_makes_a_gene_not_runnable():
    assert c2.fit("mean", _rho(-1.0, -1.0)).reason == "basal_not_positive"  # d = 0.5 + 0.5
    assert c2.fit("mean", _rho(-0.1, +0.5)).reason == "inhibitor_strength_above_1"  # e > q
    assert c2.fit("mean", _rho(+0.2, +0.1)).reason == "no_activating_rule"
    equal = {"a": 0.6, "b": 0.6}
    assert c2.fit("mean", equal).ok
    assert c2.fit("competition", equal).reason == "basal_not_positive"  # stricter than the mean
    with pytest.raises(ValueError):
        c2.fit("max", equal)


# ---- reading the gated program and running it ------------------------------------------------------------


def test_gene_rules_reads_each_rule_its_fold_states_and_compiled_cell():
    _m, genes = _genes()
    assert set(genes) == {"TGT", "TWO"}
    tgt = genes["TGT"]
    assert tgt.elements == ["E1", "E2", "E3"]
    assert tgt.rho["E1"] == pytest.approx(2**-0.5, rel=1e-3)
    assert tgt.compiled == {"E1": "placenta", "E2": "K562", "E3": "liver"}
    k562, hepg2, gm = (CELLS.index(c) for c in ("K562", "HepG2", "GM12878"))
    assert [tgt.states[e][k562] for e in tgt.elements] == ["open", "open", "closed"]
    assert [tgt.states[e][hepg2] for e in tgt.elements] == ["open", "closed", "open"]
    assert [tgt.states[e][gm] for e in tgt.elements] == ["closed", "open", "closed"]
    assert {tgt.states[e][CELLS.index("testis")] for e in tgt.elements} == {"unknown"}
    assert c2.classify(tgt)[1] == ""


@pytest.mark.parametrize("form", c2.FORMS)
def test_the_runtime_gives_the_closed_form_in_every_reader_cell_and_every_removal(form):
    m, genes = _genes()
    for gene, gr in genes.items():
        fits, why = c2.classify(gr)
        assert why == ""
        f = fits[form]
        assert f.ok, (gene, f.reason)
        out = c2.run_gene_form(rg.gated_gene_module(m, gene), gr, f, 10.0, LABEL)
        assert out["harness_miss"] == []
        assert max(out["s"]) - min(out["s"]) > 0.05  # the reader makes the cells differ
        for a, b in zip(out["s"], out["closed_s"], strict=True):
            assert a == pytest.approx(b, abs=1e-9)


def test_a_held_activator_is_counted_in_the_mean_by_the_runtime():
    """The structural fact the mean's joint fit rests on, read off the engine: in GM12878 only E2 is
    open, so the run equals basal + E2's own contribution (not E2 alone averaged over one rule)."""
    m, genes = _genes()
    gr = genes["TGT"]
    f = c2.fit("mean", gr.rho)
    small = c2.fitted_module(rg.gated_gene_module(m, "TGT"), "TGT", f, 1.0, LABEL)
    present = dict.fromkeys(gr.elements, 1.0)
    lv = c2._run_level(small, "TGT", {"cell_type": "GM12878"}, present, "mean")
    d2 = 1 - gr.rho["E2"]
    boost = 1 + (gr.rho["E3"] - 1) / f.q  # E3 (inhibiting) is closed there
    assert lv == pytest.approx((1 - f.q) + d2 * boost, rel=1e-9)
    assert small.entities["TGT"].attrs[grn.RATE_PROVENANCE_ATTRS[2]] == "t"


def test_competition_writes_the_half_saturation_the_runtime_reads():
    m, genes = _genes()
    gr = genes["TWO"]
    f = c2.fit("competition", gr.rho)
    small = c2.fitted_module(rg.gated_gene_module(m, "TWO"), "TWO", f, 1.0, LABEL)
    assert small.entities["TWO"].attrs[grn.PROMOTER_HALF_SATURATION] == pytest.approx(f.half_saturation)
    assert "competition" not in grn.ACTIVATOR_RULES and grn.ACTIVATOR_COMBINATION == "mean"


def test_multi_rule_genes_counts_rule_lines_over_programs():
    other = "\n".join(
        ["module human.noncoding.chrU", _element("H1", 100, "TWO", -0.2), _rule("H1", "TWO", -0.2, "K562")]
    )
    measured = (
        PROGRAM
        + "rule E9_measured activates ONE { strength: 0.2; when: cell_type = K562; evidence: experimental"
        ' "CRISPRi, silenced in K562"; confidence: 0.90 }\n'
    )
    out, why = c2.multi_rule_genes({"chrT": measured, "chrU": other})
    assert out == {"TGT": "chrT"}
    assert why["rules_in_several_programs"] == 1  # TWO
    assert why["a_rule_not_predicted"] == 1  # ONE, with its measured rule
    assert why["genes_with_any_rule"] == 3


def test_the_compiled_context_census_fits_only_pairs_observed_in_one_track():
    text = PROGRAM + "\n".join([_element("F3", 1300, "TWO", -0.2), _rule("F3", "TWO", -0.2, "K562")])
    out = c2.compiled_context_census({"chrT": text})
    # only TWO in K562 (F1, F3) has two mechanisms in one compiled context
    assert out["totals"]["pairs_two_or_more_mechanisms"] == 1
    assert out["totals"]["runnable_mean"] == 1
    assert out["per_context"]["K562"]["pairs_two_or_more_mechanisms"] == 1


# ---- scoring ---------------------------------------------------------------------------------------------


def _synthetic_rows(form: str, n: int, truth: bool, seed: int = 7) -> list[dict]:
    rng = random.Random(seed)
    rows = []
    while len(rows) < n:
        rho = {f"g{len(rows)}e{i}": 2.0 ** rng.uniform(-0.4, -0.05) for i in range(3)}
        rho[f"g{len(rows)}e3"] = 2.0 ** rng.uniform(0.02, 0.08)
        f = c2.fit(form, rho)
        if not f.ok:
            continue
        states = {el: [rng.choice(("open", "closed")) for _ in CELLS] for el in rho}
        s = [math.log2(c2.level(f, {el: states[el][j] == "open" for el in rho})) for j in range(len(CELLS))]
        p = rg._centre(s)
        d = [x + rng.gauss(0, 0.05) for x in p] if truth else [rng.gauss(0, 1) for _ in CELLS]
        share = [sum(states[el][j] == "open" for el in rho) / 4 for j in range(len(CELLS))]
        rows.append(
            {
                "gene": f"G{len(rows)}",
                "elements": list(rho),
                "states": states,
                "weights": f.weights,
                "e": f.e,
                "q": f.q,
                "z": f.z,
                "P": p,
                "D": rg._centre(d),
                "n_rules": 4,
                "n_inhibiting": 1,
                "open_share": rg._centre(share),
            }
        )
    return rows


@pytest.mark.parametrize("form", c2.FORMS)
def test_the_vectorised_control_is_the_scalar_closed_form(form):
    rows = _synthetic_rows(form, 30, truth=False)
    arr = c2._closed_arrays(rows, form)
    assert c2.closed_dmse(arr, arr["open"]) == pytest.approx(c2.dmse(rows), rel=1e-9)


def test_a_true_prediction_beats_and_a_random_one_does_not():
    good = c2.score_form(_synthetic_rows("mean", 120, truth=True), "mean")
    assert good["outcome"] == "BEATS_AND_CONTROL"
    assert good["reading"].startswith("the runtime's mean: beats baseline")
    bad = c2.score_form(_synthetic_rows("competition", 120, truth=False), "competition")
    assert bad["outcome"] in ("WORSE", "NO_DIFFERENCE")
    small = c2.score_form(_synthetic_rows("mean", 20, truth=True), "mean")
    assert small["outcome"] == "NOT_TESTABLE" and "fewer than the registered 100" in small["reading"]


def test_the_head_to_head_is_positive_when_competition_is_closer():
    truth = {r["gene"]: r for r in _synthetic_rows("competition", 120, truth=True)}
    noise = {g: {**r, "P": rg._centre([x * 3 for x in r["P"]])} for g, r in truth.items()}
    h = c2.head_to_head(noise, truth, 500, 1)
    assert h["compared"] and h["mse_mean_minus_mse_competition"] > 0
    assert h["clause"] == c2.HEAD_TO_HEAD["closer"]
    assert c2.head_to_head(dict(list(noise.items())[:50]), truth, 50, 1) == {"genes": 50, "compared": False}


def test_every_outcome_sentence_fills_from_the_registered_slots_only():
    nums = {"form": "F", "d": "+0.1", "lo": "+0.01", "hi": "+0.2", "n": "1,300", "g": "100", "p": "0.001"}
    for kind in ("BEATS_AND_CONTROL", "BEATS_NOT_CONTROL", "NO_DIFFERENCE", "WORSE"):
        text = c2.reading(kind, **nums)
        assert "{" not in text and "97.5% CI" in text
    assert "{" not in c2.reading("NOT_TESTABLE", form="F", n=3, min=c2.MIN_GENES)
    assert "{" not in c2.reading("HARNESS_DEFECT", bad=1, n=2)
    per = {f: {"reading": f"{f} sentence"} for f in c2.FORMS}
    out = c2.combined_reading(per, {"genes": 5, "compared": False})
    assert out.startswith("mean sentence; competition sentence. Head-to-head not compared")


def test_the_registered_words_state_the_independence_caveat_and_claim_no_function():
    assert c2.LIMITATIONS[0].startswith("independence:")
    assert "trained on" in c2.LIMITATIONS[0]
    assert any("not evidence that it does" in x for x in c2.LIMITATIONS)
    assert "Bothma" in c2.FORM_SOURCES["competition"] and "eqs. 7-9" in c2.FORM_SOURCES["competition"]
    assert c2.THETA == 0.5 and c2.FORMS == ("mean", "competition")
    assert c2.FORM_INTERVAL == (1.25, 98.75)


def test_no_outcome_sentence_strengthens():
    for s in [*c2.OUTCOMES.values(), *c2.HEAD_TO_HEAD.values()]:
        for word in ("proves", "validated", "equal", "mechanism", "function"):
            assert word not in s.lower()


def test_the_committed_registration_matches_the_code_that_applies_it():
    reg_path = Path("data/results/combine2_registration.json")
    # tracked by git from its registration commit, so its absence is a broken checkout, never a skip
    tp.must_be_present(reg_path, was="registration not yet written")
    reg = json.loads(reg_path.read_text())
    src = Path("scripts/combine2_run.py").read_text()
    bound = ast.literal_eval("{" + src.split("BOUND = {", 1)[1].split("\n}\n", 1)[0] + "\n}")
    assert len(bound) >= 30, "the run script binds too few fields"
    for field, const in bound.items():
        assert reg[field] == json.loads(json.dumps(getattr(c2, const))), field
