# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2f (lane-multrna): the multiplicative composition against measured RNA, on 7d's design.

Synthetic programs and synthetic rows only, except the registration check, which reads the committed
registration (a tracked file). No RNA value of any gene is read here.
"""

from __future__ import annotations

import math
import random

import pytest

from genomeos.attribution import combine2 as c2
from genomeos.attribution import compile as cp
from genomeos.attribution import multrna as mr
from genomeos.attribution import readergate as rg
from genomeos.attribution import saturate as sa
from genomeos.genome.reader import PeakIndex
from genomeos.lang.parser import parse

#: The tests below reach `scipy`, by one of the two doors this project opens it through:
#: `combine2.solve_z` (`scipy.optimize.brentq`, the competition form's root find) and
#: `runrules._spearman` (`scipy.stats.spearmanr`, which every registered score reads). No bare
#: environment holds it. The reasoning for naming `compose` -- the extra CI's `test` job syncs, of the
#: four that reach scipy -- and what naming one of four costs is written out at tests/test_combine2.py's
#: `needs_scipy`. The marker moves WHERE these tests run, never what they claim.
needs_scipy = pytest.mark.requires_extra("compose")

CELLS = list(mr.CELLS)

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


#: MUL: three activating elements whose single-removal drops sum past 1 (7d's fits leave no basal rate)
#: and one inhibiting element; TWO: two small activating elements every form runs; BIG: an inhibitor
#: whose fold is above 2, which the multiplicative fit refuses
ELEMENTS = {
    "M1": (100, "MUL", -1.0, "K562"),
    "M2": (300, "MUL", -0.8, "HepG2"),
    "M3": (500, "MUL", -0.6, "liver"),
    "M4": (700, "MUL", +0.3, "K562"),
    "T1": (900, "TWO", -0.3, "K562"),
    "T2": (1100, "TWO", -0.2, "HepG2"),
    "B1": (1300, "BIG", -0.4, "K562"),
    "B2": (1500, "BIG", +1.5, "K562"),
}
PROGRAM = (
    "\n".join(
        [
            "module human.noncoding.chrT",
            *(
                f'gene {g} {{ symbol: {g}; evidence: curated "GENCODE v50" symbol only, a stub }}'
                for g in ("MUL", "TWO", "BIG")
            ),
            *(_element(e, s, g, x) for e, (s, g, x, _c) in ELEMENTS.items()),
            *(_rule(e, g, x, c) for e, (_s, g, x, c) in ELEMENTS.items()),
        ]
    )
    + "\n"
)
#: biosample -> the elements a narrowPeak overlaps there; testis has no file (every element unknown)
OPEN = {
    "K562": {"M1", "M2", "M4", "T1", "B1"},
    "HepG2": {"M1", "M3", "T2", "B2"},
    "GM12878": {"M2"},
    "IMR-90": {"M3", "M4", "T1", "T2"},
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
    text, _ = cp.gate_program(PROGRAM, "chrT", selection=mr.SELECTION, readers=FakeReaders())
    m = parse(text, "t")
    genes = c2.gene_rules(m, {"MUL", "TWO", "BIG"}, "chrT", rg.compiled_cells(PROGRAM))
    return m, {gr.gene: gr for gr in genes}


# ---- the form on one gene ---------------------------------------------------------------------------------


def test_the_closed_form_is_the_product_of_the_held_elements_folds():
    _m, genes = _genes()
    gr = genes["MUL"]
    f = sa.fit(gr.rho)
    assert f.ok
    s = mr.predicted_log2(f, gr)
    for j in range(len(CELLS)):
        held = [el for el in gr.elements if gr.states[el][j] != "open"]
        assert s[j] == pytest.approx(sum(math.log2(gr.rho[el]) for el in held), abs=1e-12)
    assert max(s) - min(s) > 1.0


def test_the_runtime_gives_the_closed_form_in_every_reader_cell_and_every_removal():
    m, genes = _genes()
    for gene in ("MUL", "TWO"):
        gr = genes[gene]
        f = sa.fit(gr.rho)
        out = mr.run_gene(rg.gated_gene_module(m, gene), gr, f, 10.0, LABEL)
        assert out["harness_miss"] == []
        assert max(out["s"]) - min(out["s"]) > 0.05  # the reader makes the cells differ
        for a, b in zip(out["s"], out["closed_s"], strict=True):
            assert a == pytest.approx(b, abs=1e-9)


@needs_scipy
def test_a_gene_7d_cannot_fit_is_run_and_an_inhibitor_above_2_is_not():
    _m, genes = _genes()
    fits, why = c2.classify(genes["MUL"])
    assert why == ""
    assert sum(1 - r for r in genes["MUL"].rho.values() if r < 1) > 1
    assert {f: fits[f].reason for f in c2.FORMS} == dict.fromkeys(c2.FORMS, "basal_not_positive")
    assert sa.fit(genes["MUL"].rho).ok
    fits, why = c2.classify(genes["TWO"])
    assert why == "" and all(fits[f].ok for f in c2.FORMS) and sa.fit(genes["TWO"].rho).ok
    assert c2.classify(genes["BIG"])[1] == ""
    assert sa.fit(genes["BIG"].rho).reason == "inhibitor_strength_above_1"


def test_one_held_activator_gives_its_own_fold_under_both_and_several_fall_less_than_the_mean():
    """STRUCTURE['what_differs_from_7d'], on a gene both forms run."""
    rho = {"a": 0.7, "b": 0.8, "c": 0.9}
    mean, mult = c2.fit("mean", rho), sa.fit(rho)
    every = dict.fromkeys(rho, True)
    for el in rho:
        assert c2.level(mean, {**every, el: False}) == pytest.approx(rho[el])
        assert sa.level(mult, {**every, el: False}) == pytest.approx(rho[el])
    none = dict.fromkeys(rho, False)
    assert sa.level(mult, none) == pytest.approx(0.7 * 0.8 * 0.9)
    assert sa.level(mult, none) > c2.level(mean, none) == pytest.approx(1 - 0.6)


# ---- the score ---------------------------------------------------------------------------------------------


def _rows(n: int, truth: bool, seed: int = 7) -> list[dict]:
    rng = random.Random(seed)
    rows = []
    for g in range(n):
        k = rng.randint(2, 6)
        rho = {f"g{g}e{i}": 2.0 ** rng.uniform(-1.2, -0.05) for i in range(k)}
        rho[f"g{g}i"] = 2.0 ** rng.uniform(0.02, 0.5)
        states = {el: [rng.choice(("open", "closed")) for _ in CELLS] for el in rho}
        f = sa.fit(rho)
        s = [math.log2(sa.level(f, {el: states[el][j] == "open" for el in rho})) for j in range(len(CELLS))]
        p = rg._centre(s)
        d = [x + rng.gauss(0, 0.05) for x in p] if truth else [rng.gauss(0, 1) for _ in CELLS]
        share = [sum(states[el][j] == "open" for el in rho) / len(rho) for j in range(len(CELLS))]
        rows.append(
            {
                "gene": f"G{g}",
                "elements": list(rho),
                "states": states,
                "log2_rho": {el: math.log2(r) for el, r in rho.items()},
                "P": p,
                "D": rg._centre(d),
                "n_rules": len(rho),
                "n_inhibiting": 1,
                "open_share": rg._centre(share),
            }
        )
    return rows


def test_the_vectorised_control_is_the_scalar_closed_form():
    rows = _rows(30, truth=False)
    arr = mr._closed_arrays(rows)
    assert mr.closed_dmse(arr, arr["open"]) == pytest.approx(c2.dmse(rows), rel=1e-9)


@needs_scipy
def test_a_true_prediction_beats_and_a_random_one_does_not():
    good = mr.score_form(_rows(120, truth=True))
    assert good["outcome"] == "BEATS_AND_CONTROL"
    assert good["reading"].startswith("the multiplicative composition (Zhou et al. 2024): beats baseline")
    assert good["closed_form_delta_mse"] == pytest.approx(good["delta_mse"], rel=1e-9)
    bad = mr.score_form(_rows(120, truth=False))
    assert bad["outcome"] in ("WORSE", "NO_DIFFERENCE")
    small = mr.score_form(_rows(20, truth=True))
    assert small["outcome"] == "NOT_TESTABLE" and "fewer than the registered 100" in small["reading"]


def _by_gene(rows: list[dict]) -> dict[str, dict]:
    return {r["gene"]: r for r in rows}


def test_the_head_to_heads_are_positive_when_the_multiplicative_form_is_closer():
    truth = _by_gene(_rows(120, truth=True))
    noisy = {g: {**r, "P": rg._centre([x * 3 for x in r["P"]])} for g, r in truth.items()}
    flat = {g: {**r, "P": [0.0] * len(CELLS)} for g, r in truth.items()}
    h = mr.head_to_heads({mr.FORM: truth, "mean": noisy, "competition": flat}, 500, 1)
    assert h["compared"] and h["genes"] == 120 and h["pairs"] == 120 * len(CELLS)
    for f in c2.FORMS:
        assert h[f]["mse_form_minus_mse_multiplicative"] > 0 and h[f]["kind"] == "closer"
    assert h["mean"]["clause"] == "closer to measured RNA than the runtime's mean"
    reversed_ = mr.head_to_heads({mr.FORM: noisy, "mean": truth, "competition": truth}, 500, 1)
    assert reversed_["mean"]["kind"] == "further"
    fewer = {g: truth[g] for g in list(truth)[:50]}
    assert mr.head_to_heads({mr.FORM: fewer, "mean": truth, "competition": truth}, 50, 1) == {
        "genes": 50,
        "compared": False,
    }
    other_d = {g: {**r, "D": [x + 1 for x in r["D"]]} for g, r in truth.items()}
    with pytest.raises(ValueError):
        mr.head_to_heads({mr.FORM: truth, "mean": other_d, "competition": truth}, 50, 1)


def test_the_reading_fills_from_the_registered_slots_only():
    nums = {"form": "F", "d": "+0.1", "lo": "+0.01", "hi": "+0.2", "n": "1,300", "g": "100", "p": "0.001"}
    for kind in ("BEATS_AND_CONTROL", "BEATS_NOT_CONTROL", "NO_DIFFERENCE", "WORSE"):
        text = mr.reading(kind, **nums)
        assert "{" not in text and "97.5% CI" in text
        assert text == c2.reading(kind, **nums)  # 7d's sentences, word for word
    assert "{" not in mr.reading("COMPARATOR_DEFECT", what="x")
    assert "{" not in mr.reading("PREDICTOR_DEFECT", what="x")
    assert "{" not in mr.reading("HARNESS_DEFECT", bad=1, n=2)
    truth = _by_gene(_rows(120, truth=True))
    h = mr.head_to_heads({mr.FORM: truth, "mean": truth, "competition": truth}, 200, 1)
    out = mr.combined_reading({"reading": "M sentence"}, h)
    assert out.startswith("M sentence. Head-to-head over the 120 genes in all three populations (1,560")
    assert "{" not in out and out.count("97.5% CI") == 2
    assert "not detectably different from the runtime's mean" in out
    out = mr.combined_reading({"reading": "M"}, {"genes": 5, "compared": False})
    assert out == (
        "M. Head-to-head not compared: 5 genes are in all three populations, fewer than the registered 100."
    )


def test_assemble_applies_7ds_row_check_and_carries_the_folds():
    rows = _rows(3, truth=True)
    runs = [
        {
            **{k: r[k] for k in ("elements", "states", "log2_rho", "n_rules", "n_inhibiting", "open_share")},
            "gene": r["gene"],
            "symbol": f"S{i}",
            "chrom": "chrT",
            "s": [x + 5.0 for x in r["P"]],
        }
        for i, r in enumerate(rows)
    ]
    ids = {"S0": "ENSG0", "S1": "ENSG1"}  # S2 names no single gene_id
    tpm = {c: {"ENSG0": 3.0 + j, "ENSG1": 1.0} for j, c in enumerate(CELLS)}
    tpm[CELLS[4]].pop("ENSG1")  # ENSG1 is missing from one table
    out, why = mr.assemble(runs, tpm, ids)
    assert [r["gene"] for r in out] == ["G0"]
    assert dict(why) == {
        "symbol_not_on_exactly_one_gene_id": 1,
        "gene_id_not_on_exactly_one_row_of_every_table": 1,
    }
    assert out[0]["log2_rho"] == rows[0]["log2_rho"] and "q" not in out[0]
    assert out[0]["P"] == pytest.approx(rows[0]["P"])
    assert out[0]["D"] == pytest.approx(rg._centre([math.log2(3.0 + j + 1.0) for j in range(len(CELLS))]))


def test_the_comparator_check_names_a_population_or_value_that_moved():
    rows = {f: [{"gene": f"G{i}", "P": [0.1] * 13, "D": [0.2] * 13} for i in range(3)] for f in c2.FORMS}
    cols = ["gene", "symbol", "gene_id", "chrom", "n_rules", "n_inhibiting", "q", "P", "D"]
    table = {
        "columns": cols,
        "rows": [[f"G{i}", "", "", "", 2, 0, 0.1, [0.1] * 13, [0.2] * 13] for i in range(3)],
    }
    committed = {"genes_tables": dict.fromkeys(c2.FORMS, table), "reading": "R"}
    assert mr.comparator_problems(rows, committed, "R") == []
    assert mr.comparator_problems(rows, committed, "R'") == [
        "7d's reading, recomputed, differs from the committed reading"
    ]
    moved = {**rows, "mean": [{**rows["mean"][0], "P": [0.1 + 2e-5] * 13}, *rows["mean"][1:]]}
    assert mr.comparator_problems(moved, committed, "R")[0].startswith("mean: P or D differs")
    short = {**rows, "competition": rows["competition"][:2]}
    assert mr.comparator_problems(short, committed, "R")[0].startswith("competition: population of 2 genes")


# ---- the registered words ----------------------------------------------------------------------------------


def test_the_registered_words_keep_7ds_design_and_claim_no_function():
    assert mr.LIMITATIONS[0] == c2.LIMITATIONS[0] and mr.LIMITATIONS[0].startswith("independence:")
    assert any("not evidence that it does" in x for x in mr.LIMITATIONS)
    assert mr.FORM == "multiplicative" and mr.COMPARATOR_FORMS == ("mean", "competition")
    assert mr.FORM_INTERVAL == c2.FORM_INTERVAL == (1.25, 98.75)
    assert mr.HEAD_TO_HEAD_INTERVAL == (1.25, 98.75)
    assert (mr.BOOTSTRAP["seed"], mr.BOOTSTRAP["resamples"]) == (
        c2.BOOTSTRAP["seed"],
        c2.BOOTSTRAP["resamples"],
    )
    assert mr.PERMUTATIONS == c2.PERMUTATIONS and mr.CONTROL_ALPHA == c2.CONTROL_ALPHA
    assert (mr.MIN_GENES, mr.PSEUDOCOUNT, mr.INFORMATIVE_LOG2) == (
        c2.MIN_GENES,
        c2.PSEUDOCOUNT,
        c2.INFORMATIVE_LOG2,
    )
    assert mr.CELLS == c2.CELLS and mr.SELECTION == c2.SELECTION and mr.RUNTIME == c2.RUNTIME
    assert {k: mr.OUTCOMES[k] for k in c2.OUTCOMES} == c2.OUTCOMES
    assert "Zhou" in mr.FORM_SOURCE and mr.FIT[mr.FORM] == sa.FIT[mr.FORM]
    assert set(mr.registered_fields()) == set(mr.BOUND)


def test_no_outcome_sentence_strengthens():
    for s in [*mr.OUTCOMES.values(), *mr.HEAD_TO_HEAD.values(), mr.CENSUS_SENTENCE]:
        for word in ("proves", "validated", "equal", "mechanism", "function"):
            assert word not in s.lower()


# ---- the run script's predictor-side steps (no RNA) ------------------------------------------------------


def _run_script():
    import sys
    from pathlib import Path

    scripts = str(Path(__file__).resolve().parents[1] / "scripts")
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    import multrna_run

    return multrna_run


def _record(gene: str, ok: bool, informative: bool, miss: list[str]) -> dict:
    fr: dict = {"ok": ok, "reason": "" if ok else "inhibitor_strength_above_1"}
    if ok:
        fr["informative"] = informative
        if informative:
            fr.update(
                {"s": [0.5] * 13, "closed_s": [0.5] * 13, "harness_miss": miss, "log2_rho": {"e": -1.0}}
            )
    return {
        "gene": gene,
        "symbol": gene,
        "chrom": "chrT",
        "n_rules": 2,
        "n_inhibiting": 0,
        "excluded": "",
        "elements": ["e", "f"],
        "states": {"e": ["open"] * 13, "f": ["closed"] * 13},
        "open_share": [0.0] * 13,
        "forms": {mr.FORM: fr},
    }


def test_the_run_takes_only_informative_multiplicative_genes_and_stops_on_a_harness_miss():
    run = _run_script()
    genes = [
        _record("A", True, True, []),
        _record("B", True, False, []),
        _record("C", False, False, []),
        {"gene": "X", "excluded": "no_activating_rule", "forms": {}},
    ]
    assert [r["gene"] for r in run.mult_runs(genes)] == ["A"]
    assert run.mult_runs(genes)[0]["log2_rho"] == {"e": -1.0}
    hc = run.harness_check([], genes)
    assert (hc["failing"], hc["multiplicative_gene_runs"], hc["max_abs_run_minus_closed_form_log2"]) == (
        0,
        1,
        0.0,
    )
    genes[0]["forms"][mr.FORM]["harness_miss"] = ["cell K562: run 1, closed form 2"]
    assert run.harness_check([], genes)["failing_examples"] == [("A", mr.FORM)]


def test_the_run_descriptives_and_table_name_their_genes():
    run = _run_script()
    rows = _rows(120, truth=True)
    by_form = {mr.FORM: rows, "mean": rows[:60], "competition": rows[30:90]}
    h = mr.head_to_heads({f: _by_gene(rs) for f, rs in by_form.items()}, 50, 1)
    assert h == {"genes": 30, "compared": False}
    d = run.descriptives(by_form, {"G100", "G101", "G999"}, h)
    assert d["new_under_the_multiplicative_form"]["genes"] == 2
    assert d["on_7d_mean_population"]["genes"] == 60
    assert d["on_the_genes_in_all_three"]["genes"] == 30
    rows[0].update({"symbol": "S", "gene_id": "ENSG", "chrom": "chrT"})
    t = run._table(rows[:1])
    assert t["columns"][-2:] == ["P", "D"] and len(t["rows"][0][-1]) == 13


# ---- the committed registration ---------------------------------------------------------------------------


def test_the_committed_registration_is_the_module_it_binds():
    import json
    from pathlib import Path

    import tracked_paths as tp

    reg_path = Path("data/results/multrna_registration.json")
    # tracked by git from its registration commit, so its absence is a broken checkout, never a skip
    tp.must_be_present(reg_path, was="registration not yet written")
    reg = json.loads(reg_path.read_text())
    assert {f: reg[f] for f in mr.BOUND} == mr.registered_fields()
    assert reg["census_facts"]["seven_d_census_reproduced"]["problems"] == []
    m = reg["census_facts"]["reader_selection"]["forms"][mr.FORM]
    assert m["run"] == m["informative"] and m["harness_failing_genes"] == 0
    assert reg["seven_d_committed_verbatim"]["outcomes"] == {"mean": "WORSE", "competition": "WORSE"}
    assert reg["result_manifest"]["complete"] is True
    assert reg["result_manifest"]["code_cleanliness"]["own_code_is_committed"] is True
    assert reg["alphagenome_requests"] == 0
    assert any(p.startswith("data/cache/readergate/rna") for p in reg["pinned_sha256"])
