# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-capsumrna: the capped sum (grn 'sum_capped') on 7d's design, the fit, its two roots and the census.

Synthetic programs only. No RNA value of any gene is read here.
"""

from __future__ import annotations

import itertools
import math
import random

import pytest

from genomeos.attribution import capsumrna as cs
from genomeos.attribution import combine2 as c2
from genomeos.attribution import compile as cp
from genomeos.attribution import readergate as rg
from genomeos.genome.reader import PeakIndex
from genomeos.lang.parser import parse
from genomeos.runtime import grn

CELLS = list(cs.CELLS)

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


#: TWO: two activating elements (W below the cap, one root); FOUR: a dominant activator, three small ones
#: and an inhibitor (two roots); PAST: three equal activators whose W at the declared strength is 1.5
ELEMENTS = {
    "T1": (100, "TWO", -0.300, "K562"),
    "T2": (300, "TWO", -0.200, "HepG2"),
    "F1": (500, "FOUR", -0.322, "K562"),
    "F2": (700, "FOUR", -0.074, "HepG2"),
    "F3": (900, "FOUR", -0.074, "liver"),
    "F4": (1100, "FOUR", -0.074, "K562"),
    "F5": (1300, "FOUR", +0.300, "K562"),
    "P1": (1500, "PAST", -0.322, "K562"),
    "P2": (1700, "PAST", -0.322, "HepG2"),
    "P3": (1900, "PAST", -0.322, "liver"),
}
GENES = ("TWO", "FOUR", "PAST")
PROGRAM = (
    "\n".join(
        [
            "module human.noncoding.chrT",
            *(
                f'gene {g} {{ symbol: {g}; evidence: curated "GENCODE v50" symbol only, a stub }}'
                for g in GENES
            ),
            *(_element(e, s, g, x) for e, (s, g, x, _c) in ELEMENTS.items()),
            *(_rule(e, g, x, c) for e, (_s, g, x, c) in ELEMENTS.items()),
        ]
    )
    + "\n"
)
#: biosample -> the elements a narrowPeak overlaps there; testis has no file (every element unknown)
OPEN = {
    "K562": {"T1", "F1", "F2", "F5", "P1", "P2"},
    "HepG2": {"T2", "F2", "F3", "F4", "P3"},
    "GM12878": {"T1", "T2", "F1", "F4", "P1"},
    "IMR-90": {"F3", "F5", "P2", "P3"},
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


def _genes():
    text, _ = cp.gate_program(PROGRAM, "chrT", selection=cs.SELECTION, readers=FakeReaders())
    m = parse(text, "t")
    genes = c2.gene_rules(m, set(GENES), "chrT", rg.compiled_cells(PROGRAM))
    return m, {gr.gene: gr for gr in genes}


def _subsets(keys):
    for bits in itertools.product((True, False), repeat=len(keys)):
        yield dict(zip(keys, bits, strict=True))


# ---- the form is the engine's, unchanged -------------------------------------------------------------------


def test_the_form_is_the_engines_sum_capped_and_the_default_stays_the_mean():
    assert cs.FORM in grn.ACTIVATOR_RULES and cs.FORM not in grn.OPT_IN_ACTIVATOR_RULES
    assert grn.ACTIVATOR_COMBINATION == "mean"
    assert grn.combine_activators([0.5, 0.7], cs.FORM) == 1.0
    assert grn.combine_activators([0.5, 0.2, 0.0], cs.FORM) == pytest.approx(0.7)
    assert grn.combine_activators([0.5, 0.2, 0.0]) == pytest.approx(0.7 / 3)
    assert cs.H == 0.5  # the bridge's H_INTACT, which with s = 1 puts the strongest term at half the cap


# ---- the fit and its two roots -----------------------------------------------------------------------------


def test_the_uncapped_root_is_the_mean_in_every_open_set():
    """STRUCTURE['identity_with_the_mean'] on random genes below the cap."""
    rnd, seen = random.Random(5), 0
    while seen < 200:
        k = rnd.randint(2, 7)
        rho = {f"e{i}": 1.0 - rnd.uniform(0.01, 0.3) for i in range(k)}
        if rnd.random() < 0.5:
            rho["inh"] = 1.0 + rnd.uniform(0.01, 0.3)
        f = cs.fit(rho)
        if not f.ok:
            assert f.reason in ("sum_past_the_cap", "basal_not_positive", "inhibitor_strength_above_1")
            continue
        seen += 1
        mean = c2.fit("mean", rho)
        assert f.w <= 1.0 + 1e-9 and f.q == mean.q and f.s == mean.s
        for o in _subsets(list(rho)):
            assert cs.level(f, o) == pytest.approx(c2.level(mean, o), abs=1e-12)


def test_past_the_cap_there_is_no_fit_on_either_root():
    """STRUCTURE['no_fit_past_the_cap']."""
    rnd, seen = random.Random(11), 0
    while seen < 300:
        k = rnd.randint(3, 9)
        rho = {f"e{i}": 1.0 - rnd.uniform(0.01, 0.2) for i in range(k)}
        if cs.H * sum(1 - r for r in rho.values()) / max(1 - r for r in rho.values()) <= 1.0:
            continue
        seen += 1
        assert cs.fit(rho).reason in ("sum_past_the_cap", "basal_not_positive")
        assert not cs.capped_root(rho).ok


def test_the_capped_root_reproduces_every_removal_and_differs_from_the_mean():
    """STRUCTURE['a_second_root']: four activators with a dominant one, and an inhibitor."""
    rho = {"a": 0.8, "b": 0.95, "c": 0.95, "d": 0.95, "i": 1.2}
    un, cap = cs.fit(rho), cs.capped_root(rho)
    assert un.ok and cap.ok and un.root == "uncapped" and cap.root == "capped"
    assert cap.w > 1.0 and un.w < 1.0
    assert max(cap.s[el] for el in "abcd") == pytest.approx(1.0)
    every = dict.fromkeys(rho, True)
    for f in (un, cap):
        assert cs.level(f, every) == pytest.approx(1.0)
        for el in rho:
            assert cs.level(f, {**every, el: False}) == pytest.approx(rho[el])
    mean = c2.fit("mean", rho)
    gaps = [abs(cs.level(cap, o) - c2.level(mean, o)) for o in _subsets(list(rho))]
    assert max(gaps) > 0.01
    assert cs.capped_root({"a": 0.8, "b": 0.95, "c": 0.95}).reason == "no_capped_root"  # k = 3


def test_classify_names_each_gene_status():
    _m, genes = _genes()
    assert c2.classify(genes["PAST"])[0]["mean"].ok  # the mean runs it; the capped sum cannot
    got = {g: cs.classify(gr)["status"] for g, gr in genes.items()}
    assert got == {"TWO": "runnable", "FOUR": "fit_not_unique", "PAST": "sum_past_the_cap"}


# ---- the runtime -------------------------------------------------------------------------------------------


def test_the_runtime_gives_the_closed_form_on_both_roots_and_the_mean_on_the_uncapped_one():
    m, genes = _genes()
    for gene in ("TWO", "FOUR"):
        gr = genes[gene]
        small = rg.gated_gene_module(m, gene)
        roots = cs.classify(gr)
        mean_run = c2.run_gene_form(small, gr, c2.fit("mean", gr.rho), 10.0, LABEL)
        for root in ("uncapped", "capped"):
            f = roots[root]
            if not f.ok:
                continue
            out = cs.run_gene(small, gr, f, 10.0, LABEL)
            assert out["harness_miss"] == []
            assert max(out["s"]) - min(out["s"]) > 0.05  # the reader makes the cells differ
            for a, b in zip(out["s"], out["closed_s"], strict=True):
                assert a == pytest.approx(b, abs=1e-9)
            same = max(abs(a - b) for a, b in zip(out["s"], mean_run["s"], strict=True))
            if root == "uncapped":
                assert same < 1e-9
                assert cs.difference_from_mean(f, gr) < 1e-12
            else:
                assert same > 1e-3
    assert not cs.classify(genes["TWO"])["capped"].ok


# ---- the decision ------------------------------------------------------------------------------------------


def _summary(diff: float | None, miss: int = 0) -> dict:
    run = {"run": 3, "informative": 3, "harness_failing_genes": miss}
    return {
        "unique_fit": {"max_abs_difference_from_mean_log2": diff},
        "roots": {"uncapped": run, "capped": run},
    }


def test_the_decision_registers_only_when_a_unique_fit_differs_from_the_mean_or_the_runtime_misses():
    assert cs.decide(_summary(1e-15)) == "no_registration_identity_with_the_mean"
    assert cs.decide(_summary(None)) == "no_registration_identity_with_the_mean"
    assert cs.decide(_summary(0.2)) == "register_on_the_genes_that_differ"
    assert cs.decide(_summary(1e-15, miss=1)) == "register_on_the_genes_that_differ"


def test_no_word_claims_function_or_strengthens():
    words = " ".join([cs.CENSUS, cs.CENSUS_SENTENCE, cs.DECISION, *cs.STRUCTURE.values(), *cs.FIT.values()])
    for bad in ("explains expression", "explain expression", "regulates", "causal", "proves"):
        assert bad not in words
    assert "7d's run has already read the RNA" in cs.DECISION
    assert math.isclose(cs.CLOSED_FORM_TOLERANCE, c2.CLOSED_FORM_TOLERANCE)
