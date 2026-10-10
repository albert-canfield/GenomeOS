# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1 (lane-runrules): the compiled rules, run, against measured RNA. Fixtures only; no cached data read.

What is pinned: the panel mapping, that a run of a bridged rule returns its input observation read as
a present/absent ratio (and nothing that depends on the rate), which pairs are admitted and why the
rest are not, and that every registered outcome is reachable and worded as registered.
"""

from __future__ import annotations

import ast
import gzip
import json
import random
from pathlib import Path

import pytest

from genomeos.attribution import bridge
from genomeos.attribution import runrules as rr
from genomeos.lang.parser import parse


def _rule(el: str, action: str, gene: str, cell: str, x: float) -> str:
    """A rule as the compiler writes a predicted one (compile.compile_chromosome)."""
    return (
        f"rule {el} {action} {gene} {{ strength: {abs(x)}; when: cell_type = {cell}; "
        f'evidence: predicted "AlphaGenome deletion, {cell}" effect {x:+.1f} log2 fold change on deletion,'
        " probability unavailable }\n"
    )


PROGRAM = (
    "gene GA { symbol: GA }\ngene GB { symbol: GB }\ngene GC { symbol: GC }\n"
    + "".join(f"element E{i} {{ class: enhancer }}\n" for i in range(1, 5))
    + _rule("E1", "activates", "GA", "Whole_Blood", -1.0)
    + _rule("E2", "inhibits", "GB", "Liver", 0.5)
    + _rule("E3", "activates", "GC", "Lung", -0.3)
    + _rule("E4", "activates", "GC", "Lung", -0.2)
)
BORROWED = 1.76


@pytest.fixture
def program():
    return parse(PROGRAM, "fixture")


# ---- the panel ---------------------------------------------------------------------------------------


def test_panel_is_the_54_gtex_tissues_by_exact_name():
    assert len(rr.GTEX_V8_TISSUES) == 54
    assert len(rr.PANEL) == 54  # no two GTEx names collapse to one label
    assert rr.tissue_label("Brain - Spinal cord (cervical c-1)") == "Brain_Spinal_cord_cervical_c_1"
    assert rr.tissue_label("Cells - EBV-transformed lymphocytes") == "Cells_EBV_transformed_lymphocytes"
    assert rr.tissue_label("Skin - Not Sun Exposed (Suprapubic)") == "Skin_Not_Sun_Exposed_Suprapubic"
    assert rr.tissue_label("Whole Blood") == "Whole_Blood"
    # case-sensitive: the ENCODE biosample `liver` is not the GTEx track `Liver`
    assert "Liver" in rr.PANEL and "liver" not in rr.PANEL


# ---- the run returns its input, read as present against absent ---------------------------------------


def test_a_run_of_an_activating_rule_gives_minus_log2_rho(program):
    pair = rr.run_pair(rr.one_gene_module(program, "GA"), "GA", {"cell_type": "Whole_Blood"}, {}, BORROWED)
    assert pair is not None
    assert pair["action"] == "activates"
    assert pair["rho"] == pytest.approx(0.5)
    assert pair["m"] == pytest.approx(1.0, abs=rr.CLOSED_FORM_TOLERANCE)
    assert abs(pair["m"] - pair["closed_form"]) <= rr.CLOSED_FORM_TOLERANCE
    assert pair["tier"] == "borrowed_median"
    # STRUCTURE['level']: the intact level is T / delta_m with the default half-life (delta_m = 1/h)
    assert pair["L_present"] == pytest.approx(BORROWED, rel=1e-9)
    # STRUCTURE['on_off']: the element-absent level is still above zero
    assert pair["L_absent"] > 0


def test_a_run_of_an_inhibiting_rule_is_negative(program):
    pair = rr.run_pair(rr.one_gene_module(program, "GB"), "GB", {"cell_type": "Liver"}, {}, BORROWED)
    assert pair is not None and pair["action"] == "inhibits"
    assert pair["m"] == pytest.approx(-0.5, abs=rr.CLOSED_FORM_TOLERANCE)


def test_the_run_does_not_depend_on_the_rate(program):
    ctx = {"cell_type": "Whole_Blood"}
    small = rr.one_gene_module(program, "GA")
    lent = rr.run_pair(small, "GA", ctx, {}, BORROWED)
    own = rr.run_pair(small, "GA", ctx, {"GA": 250.0}, BORROWED)
    assert own["tier"] == "measured" and lent["tier"] == "borrowed_median"
    assert own["L_present"] == pytest.approx(250.0, rel=1e-9)
    assert own["m"] == pytest.approx(lent["m"], abs=1e-12)


def test_several_mechanisms_are_not_run(program):
    assert rr.run_pair(rr.one_gene_module(program, "GC"), "GC", {"cell_type": "Lung"}, {}, BORROWED) is None


def test_a_rule_outside_its_cell_is_not_run(program):
    assert rr.run_pair(rr.one_gene_module(program, "GA"), "GA", {"cell_type": "Liver"}, {}, BORROWED) is None


def test_one_gene_module_fits_what_the_whole_program_fits(program):
    for gene, cell in (("GA", "Whole_Blood"), ("GB", "Liver")):
        ctx = {"cell_type": cell}
        whole = bridge.parameterize(program, ctx, build=False, rates={}, borrowed_rate=BORROWED)
        small = bridge.parameterize(
            rr.one_gene_module(program, gene), ctx, build=False, rates={}, borrowed_rate=BORROWED
        )
        assert {k: v for k, v in whole.strengths.items() if k[1] == gene} == small.strengths


# ---- census -----------------------------------------------------------------------------------------


def test_census_counts_per_context(program):
    c = rr.census_module(program, {}, BORROWED, {}, {})
    assert c["rules"] == 4 and c["rules_with_gene_source"] == 0
    wb, lung = c["contexts"]["Whole_Blood"], c["contexts"]["Lung"]
    assert wb["single_mechanism_pairs"] == 1 and wb["resolved_borrowed"] == 1 and wb["resolved_measured"] == 0
    assert wb["unresolved_measured"] == {"gene_rate_unmeasured": 1}
    assert wb["resolved_borrowed_basal_not_positive"] == 0
    assert lung["several_mechanism_pairs"] == 1 and lung["resolved_borrowed"] == 0
    assert lung["unresolved_borrowed"] == {"not_identifiable": 1}
    pooled = rr.pool_census({"a": c, "b": c})
    assert pooled["totals"]["rules"] == 8
    assert rr.summarise(pooled["contexts"], {"Lung"})["several_mechanism_pairs"] == 2


def test_a_gene_source_is_counted():
    m = parse(
        "gene GA { symbol: GA }\ngene GB { symbol: GB }\n"
        "rule GA activates GB { strength: 0.5; when: cell_type = Liver }\n",
        "g",
    )
    assert rr.census_module(m, {}, BORROWED, {}, {})["rules_with_gene_source"] == 1


def test_rule_tissues_are_panel_contexts_only(program):
    t = rr.rule_tissues(program)
    assert t["GA"] == {"Whole_Blood"} and t["GC"] == {"Lung"}
    m = parse(
        "element E { class: enhancer }\ngene G { symbol: G }\n"
        "rule E activates G { strength: 0.5; when: cell_type = K562 }\n",
        "k",
    )
    assert rr.rule_tissues(m) == {}


# ---- measured RNA, population and score -------------------------------------------------------------


def _gct(path: Path, header: list[str], rows: list[tuple[str, list[float]]]) -> Path:
    with gzip.open(path, "wt") as fh:
        fh.write("#1.2\n")
        fh.write(f"{len(rows)}\t{len(header)}\n")
        fh.write("\t".join(["Name", "Description", *header]) + "\n")
        for i, (sym, vals) in enumerate(rows):
            fh.write("\t".join([f"ENSG{i:011d}", sym, *map(str, vals)]) + "\n")
    return path


def test_load_gtex_refuses_another_header(tmp_path):
    p = _gct(tmp_path / "x.gct.gz", ["Liver", "Lung"], [("GA", [1.0, 2.0])])
    with pytest.raises(ValueError, match="not the registered one"):
        rr.load_gtex(p)


def test_load_gtex_drops_ambiguous_symbols(tmp_path):
    n = len(rr.GTEX_V8_TISSUES)
    rows = [("GA", [1.0] * n), ("GB", [2.0] * n), ("GB", [3.0] * n)]
    tissues, tpm, counts = rr.load_gtex(_gct(tmp_path / "x.gct.gz", list(rr.GTEX_V8_TISSUES), rows))
    assert tissues == list(rr.GTEX_V8_TISSUES)
    assert set(tpm) == {"GA"} and counts["symbols_on_several_rows"] == 1


def test_assemble_computes_D_and_counts_exclusions():
    tissues = list(rr.GTEX_V8_TISSUES)
    col = {rr.tissue_label(t): i for i, t in enumerate(tissues)}
    row = [1.0] * 54
    row[col["Liver"]] = 15.0  # log2(16) - log2(2) = 3
    tpm = {"GA": row, "GB": [0.0] * 54}
    pairs = [
        {"gene": "GA", "tissue": "Liver", "m": 0.5},
        {"gene": "GB", "tissue": "Liver", "m": 0.5},
        {"gene": "GX", "tissue": "Liver", "m": 0.5},
    ]
    many = {"GB": set(list(col)[:40])}  # 14 free tissues, fewer than 20
    pop, why = rr.assemble(pairs, {"GA": {"Liver"}, **many}, tpm, tissues)
    assert [p["gene"] for p in pop] == ["GA"]
    assert pop[0]["D"] == pytest.approx(3.0)
    assert pop[0]["on_in_rule_tissue"] is True and pop[0]["on_share_reference"] == 1.0
    assert why == {"symbol_not_on_exactly_one_gtex_row": 1, "fewer_than_min_reference_tissues": 1}


def _pop(make_d, n_genes=400, tissues=("Liver", "Lung", "Whole_Blood"), seed=1):
    rng = random.Random(seed)
    pop = []
    for g in range(n_genes):
        for t in tissues[: 1 + g % 2]:
            m = rng.choice([1, 1, 1, -1]) * rng.uniform(0.05, 1.0)
            pop.append(
                {
                    "gene": f"G{g}",
                    "tissue": t,
                    "m": m,
                    "action": "activates" if m > 0 else "inhibits",
                    "tier": "measured" if g % 3 == 0 else "borrowed_median",
                    "on_in_rule_tissue": True,
                    "on_share_reference": 0.5,
                    "D": make_d(m, rng),
                }
            )
    return pop


def test_score_reaches_each_registered_outcome():
    assert rr.score(_pop(lambda m, r: 2 * m + r.gauss(0, 0.3)))["outcome"] == "BEATS_AND_CONTROL"
    # a shift every gene of a tissue shares (the gate), not tied to the gene's own rule
    assert rr.score(_pop(lambda m, r: 2.0 + r.gauss(0, 0.3)))["outcome"] == "BEATS_NOT_CONTROL"
    # an uninformative shift costs m^2, so the null that balances is D = m / 2 plus noise
    assert rr.score(_pop(lambda m, r: 0.5 * m + r.gauss(0, 3.0)))["outcome"] == "NO_DIFFERENCE"
    assert rr.score(_pop(lambda m, r: r.gauss(0, 0.3)))["outcome"] == "WORSE"
    assert rr.score(_pop(lambda m, r: -2 * m + r.gauss(0, 0.3)))["outcome"] == "WORSE"


def test_score_stops_below_the_registered_size():
    out = rr.score(_pop(lambda m, r: m, n_genes=100))
    assert out["outcome"] == "NOT_TESTABLE" and out["pairs"] < rr.MIN_POPULATION
    assert out["reading"].startswith("not testable at the registered size")


def test_score_descriptives_add_up():
    out = rr.score(_pop(lambda m, r: 2 * m + r.gauss(0, 0.3)))
    d = out["descriptives"]
    assert sum(d["sign_table"].values()) == out["pairs"]
    assert d["gate_only_agreement_D_positive"] + d["runtime_literal_agreement_D_negative"] <= 1.0
    assert sum(v["pairs"] for v in d["by_tissue"].values()) == out["pairs"]
    assert out["mse_baseline"] - out["mse_run"] == pytest.approx(out["delta_mse"])
    assert out["limitations"] == list(rr.LIMITATIONS)


def test_delta_mse_is_the_registered_formula():
    m, d = [0.5, -0.2], [1.0, 0.3]
    base = sum(x * x for x in d) / 2
    run = sum((x - y) ** 2 for x, y in zip(d, m, strict=True)) / 2
    assert rr.delta_mse(m, d) == pytest.approx(base - run)


def test_decide_boundaries():
    assert rr.decide(0.1, 0.01, 0.2, 0.01) == "BEATS_AND_CONTROL"
    assert rr.decide(0.1, 0.01, 0.2, 0.5) == "BEATS_NOT_CONTROL"
    assert rr.decide(0.0, -0.01, 0.01, 0.5) == "NO_DIFFERENCE"
    assert rr.decide(0.0, 0.0, 0.01, 0.01) == "NO_DIFFERENCE"  # a bound AT zero is not above it
    assert rr.decide(-0.1, -0.2, -0.01, 1.0) == "WORSE"


def test_every_outcome_sentence_fills():
    nums = {"d": "+0.1", "lo": "+0.01", "hi": "+0.2", "n": "1,000", "g": "900", "t": 40, "p": "0.01"}
    for kind in ("BEATS_AND_CONTROL", "BEATS_NOT_CONTROL", "NO_DIFFERENCE", "WORSE"):
        assert "{" not in rr.reading(kind, **nums)
    assert "{" not in rr.reading("NOT_TESTABLE", n=10, min=500)
    assert "{" not in rr.reading("HARNESS_DEFECT", tol=1e-6, bad=1, n=10)


def test_no_outcome_sentence_strengthens():
    for s in rr.OUTCOMES.values():
        for word in ("proves", "validated", "equal", "silencer", "mechanism"):
            assert word not in s.lower()


def test_the_committed_registration_matches_the_code_that_applies_it():
    reg_path = Path("data/results/runrules_registration.json")
    if not reg_path.exists():
        pytest.skip("registration not yet written")
    reg = json.loads(reg_path.read_text())
    src = Path("scripts/runrules_run.py").read_text()
    bound = ast.literal_eval("{" + src.split("BOUND = {", 1)[1].split("\n}\n", 1)[0] + "\n}")
    assert bound, "the run script binds no field"
    for field, const in bound.items():
        assert reg[field] == json.loads(json.dumps(getattr(rr, const))), field


def test_the_gtex_row_is_read_by_symbol_not_by_program_identifier():
    m = parse(
        "gene NKX2_1 { symbol: NKX2-1 }\nelement E1 { class: enhancer }\n"
        + _rule("E1", "activates", "NKX2_1", "Lung", -0.5),
        "s",
    )
    pair = rr.run_pair(rr.one_gene_module(m, "NKX2_1"), "NKX2_1", {"cell_type": "Lung"}, {}, BORROWED)
    assert pair is not None and pair["gene"] == "NKX2_1" and pair["symbol"] == "NKX2-1"
    tissues = list(rr.GTEX_V8_TISSUES)
    pop, why = rr.assemble([pair], {"NKX2_1": {"Lung"}}, {"NKX2-1": [1.0] * 54}, tissues)
    assert len(pop) == 1 and not why
    pop, why = rr.assemble([pair], {"NKX2_1": {"Lung"}}, {"NKX2_1": [1.0] * 54}, tissues)
    assert pop == [] and why == {"symbol_not_on_exactly_one_gtex_row": 1}
